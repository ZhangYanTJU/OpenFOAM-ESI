/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "basicSch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRDragModels
{
    defineTypeNameAndDebug(basicSch, 0);
    addToRunTimeSelectionTable(PDRDragModel, basicSch, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRDragModels::basicSch::basicSch
(
    const dictionary& PDRProperties,
    const compressible::RASModel& turbulence,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    PDRDragModel(PDRProperties,turbulence, rho, U, phi),
    Csu("Csu", dimless, PDRDragModelCoeffs_),
    Csk("Csk", dimless, PDRDragModelCoeffs_),
    Aw_
    (
        IOobject
        (
            "Aw",
            U_.mesh().facesInstance(),
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    ),

    CR_
    (
        IOobject
        (
            "CR",
            U_.mesh().facesInstance(),
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    ),
    nrCoef_(PDRDragModelCoeffs_.get<scalar>("nrCoef")),
    nrExp2_(PDRDragModelCoeffs_.get<scalar>("nrExp2")),
    lCoef_(PDRDragModelCoeffs_.get<scalar>("lCoef")),
    maxSchFac_(PDRDragModelCoeffs_.get<scalar>("maxSchFac")),
    subGridSchelkin_(PDRDragModelCoeffs_.get<bool>("subGridSchelkin"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::PDRDragModels::basicSch::~basicSch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::PDRDragModels::basicSch::Dcu() const
{
    tmp<volSymmTensorField> tDragDcu
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tDragDcu",
                U_.mesh().time().constant(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedSymmTensor(dimMass/dimTime/dimVolume, Zero)

        )
    );

    volSymmTensorField& DragDcu = tDragDcu.ref();

    if (on_)
    {
        const volScalarField& betav =
            U_.db().lookupObject<volScalarField>("betav");

        DragDcu =
            (0.5*rho_)*CR_*mag(U_) + (Csu*I)*betav*turbulence_.muEff()*sqr(Aw_);
    }

    return tDragDcu;
}


Foam::tmp<Foam::volScalarField> Foam::PDRDragModels::basicSch::Gk() const
{
    tmp<volScalarField> tGk
    (
        new volScalarField
        (
            IOobject
            (
                "tGk",
                U_.mesh().time().constant(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)

        )
    );

    volScalarField& Gk = tGk.ref();

    if (on_)
    {
        const volScalarField& betav =
            U_.db().lookupObject<volScalarField>("betav");

        const volSymmTensorField& CT =
            U_.db().lookupObject<volSymmTensorField>("CT");

        Gk =
            (0.5*rho_)*mag(U_)*(U_ & CT & U_)
          + Csk*betav*turbulence_.muEff()*sqr(Aw_)*magSqr(U_);

        if (subGridSchelkin_)
        {
            Gk *= schFac();
        }
    }


    return tGk;
}

Foam::tmp<Foam::volScalarField> Foam::PDRDragModels::basicSch::schFac() const
{
    const volScalarField& Su_ = U_.db().lookupObject<volScalarField>("Su");
    const volScalarField& rhou_ = U_.db().lookupObject<volScalarField>("rhou");
    const volScalarField& muu_ = U_.db().lookupObject<volScalarField>("muu");

    tmp<volScalarField> tfac
    (
        new volScalarField
        (
            IOobject
            (
                "tfac",
                U_.mesh().time().constant(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar(dimless, Zero)

        )
    );

    volScalarField& schFac = tfac.ref();

    const volScalarField& k = turbulence_.k();
    const volScalarField& epsilon = turbulence_.epsilon();

    const volScalarField up(sqrt((2.0/3.0)*k));
    const volScalarField l(lCoef_*sqrt(3.0/2.0)*up*k/epsilon);

    volScalarField Rs(Su_*l*rhou_/muu_);

    if (subGridSchelkin_)
    {
        schFac = max
        (
            1.0,
            min
            (
                maxSchFac_,
                pow(Rs, 2.0 * SchelkinExponent(nrCoef_, nrExp2_, Su_))
            )
        );
    }

    return tfac;
}


Foam::tmp<Foam::volScalarField> Foam::PDRDragModels::basicSch::SchelkinExponent
(
    const scalar nrCoef,
    const scalar nrExp,
    const volScalarField& Su
) const
{
    const fvMesh& mesh = Su.mesh();

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    const volScalarField& Nv = mesh.lookupObject<volScalarField>("Nv");
    const volSymmTensorField& nsv =
        mesh.lookupObject<volSymmTensorField>("nsv");

    tmp<volScalarField> tN
    (
        new volScalarField
        (
            IOobject
            (
                "tN",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(Nv.dimensions(), Zero)

        )
    );

    volScalarField& N = tN.ref();

    N.primitiveFieldRef() = Nv.primitiveField()*pow(mesh.V(), 2.0/3.0);

    tmp<volSymmTensorField> tns
    (
        new volSymmTensorField
        (
            IOobject
            (
                "tns",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(nsv.dimensions(), Zero)
        )
    );

    volSymmTensorField& ns = tns.ref();

    ns.primitiveFieldRef() = nsv.primitiveField()*pow(mesh.V(), 2.0/3.0);

    const volVectorField Uhat
    (
        U/(mag(U) + dimensionedScalar("Usmall", U.dimensions(), 1e-4))
    );

    const volScalarField nr(sqrt(max(N - (Uhat & ns & Uhat), scalar(1.0))));

    //Re use tN
    N.primitiveFieldRef() =
        nrCoef*((scalar(1.0) - pow(nrExp, nr))/(1.0 - nrExp) - scalar(1.0));

    return tN;

}


bool Foam::PDRDragModels::basicSch::read(const dictionary& PDRProperties)
{
    PDRDragModel::read(PDRProperties);

    PDRDragModelCoeffs_.readEntry("Csu", Csu.value());
    PDRDragModelCoeffs_.readEntry("Csk", Csk.value());

    return true;
}


void Foam::PDRDragModels::basicSch::writeFields() const
{
    Aw_.write();
    CR_.write();
}

// ************************************************************************* //
