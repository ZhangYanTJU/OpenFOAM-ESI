/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "contactHeatFluxSource.H"
#include "faMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "famSup.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(contactHeatFluxSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        contactHeatFluxSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::contactHeatFluxSource::contactHeatFluxSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const fvPatch& patch
)
:
    faceSetOption(sourceName, modelType, dict, mesh, patch),
    temperatureCoupledBase(patch, dict),
    TName_(dict.get<word>("T")),
    TprimaryName_(dict.get<word>("Tprimary")),
    Tp_(mesh.lookupObject<volScalarField>(TprimaryName_)),
    Tw1_
    (
        IOobject
        (
            "Tw1_" + sourceName,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero),
        zeroGradientFaPatchScalarField::typeName
    ),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    curTimeIndex_(-1)
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained
    fieldNames_.setSize(1, TName_);

    applied_.setSize(fieldNames_.size(), false);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fa::contactHeatFluxSource::~contactHeatFluxSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::areaScalarField> Foam::fa::contactHeatFluxSource::htc() const
{
    IOobject io
    (
        "thtc",
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<areaScalarField> thtc
    (
        new areaScalarField
        (
            io,
            regionMesh(),
            dimensionedScalar(dimPower/dimArea/dimTemperature, Zero)
        )
    );

    areaScalarField& htc = thtc.ref();

    const volScalarField::Boundary& vfb = Tp_.boundaryField();

    htc.primitiveFieldRef() =
        temperatureCoupledBase::kappa
        (
            vsm().mapInternalToSurface<scalar>(vfb)()
        )*patch().deltaCoeffs();

    if (contactRes_ != 0)
    {
        tmp<areaScalarField> tcontact
        (
            new areaScalarField
            (
                io,
                regionMesh(),
                dimensionedScalar
                (
                    "contact",
                    dimPower/dimArea/dimTemperature,
                    contactRes_
                )
            )
        );
        areaScalarField& contact = tcontact.ref();
        htc.primitiveFieldRef()  += contact.primitiveField();
    }

    return thtc;
}


void Foam::fa::contactHeatFluxSource::addSup
(
    const areaScalarField& h,
    const areaScalarField& rhoCph,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo<< name() << ": applying source to " << eqn.psi().name() << endl;

        if (curTimeIndex_ != mesh().time().timeIndex())
        {
            const volScalarField::Boundary& vfb = Tp_.boundaryField();

            Tw1_.primitiveFieldRef() =
                this->vsm().mapInternalToSurface<scalar>(vfb);

            tmp<areaScalarField> htcw = htc();

            eqn += -fam::Sp(htcw(), eqn.psi()) + htcw()*Tw1_;

            curTimeIndex_ = mesh().time().timeIndex();
        }
    }
}


bool Foam::fa::contactHeatFluxSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("T", TName_);

        if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
        {
            dict.readEntry("kappaLayers", kappaLayers_);

            if (thicknessLayers_.size() > 0)
            {
                // Calculate effective thermal resistance by harmonic averaging
                forAll(thicknessLayers_, iLayer)
                {
                    contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
                contactRes_ = 1.0/contactRes_;
            }
        }

        return true;
    }

    return false;
}

// ************************************************************************* //
