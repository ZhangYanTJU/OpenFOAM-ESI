/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "buoyancyTurbSource.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(buoyancyTurbSource, 0);
    addToRunTimeSelectionTable(option, buoyancyTurbSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::calcB()
{
    const auto& alphat = mesh_.lookupObjectRef<volScalarField>("alphat");
    const auto& T = mesh_.lookupObjectRef<volScalarField>("T");
    const auto& g = meshObjects::gravity::New(mesh_.time());

    B_ = -beta_*alphat()*(fvc::grad(T) & g)();
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceOmega
(
    fvMatrix<scalar>& eqn
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    tmp<volScalarField> tomega = turbPtr->omega();
    if (tomega.isTmp())
    {
        FatalErrorInFunction
            << "Unable to find an omega field." << nl
            << "buoyancyTurbSource needs an omega field."
            << exit(FatalError);
    }

    const volScalarField::Internal& nut = turbPtr->nut()();
    const auto& gamma =
        mesh_.lookupObjectRef<volScalarField::Internal>
        (
            word(turbPtr->type() + ":gamma")
        );

    eqn += gamma/(nut + dimensionedScalar(nut.dimensions(), SMALL))*B_;
}


void Foam::fv::buoyancyTurbSource::buoyancyTurbSourceK
(
    fvMatrix<scalar>& eqn
) const
{
    eqn += B_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::buoyancyTurbSource::buoyancyTurbSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    beta_
    (
        dimensionedScalar
        (
            dimless/dimTemperature,
            coeffs_.getCheckOrDefault<scalar>
            (
                "beta",
                3.3e-3,
                [&](const scalar x){ return x > SMALL; }
            )
        )
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero)
    )
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.setSize(2, "undefined");

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (!tepsilon.isTmp())
    {
        fieldNames_[0] = tepsilon().name();
    }
    else if (!tomega.isTmp())
    {
        fieldNames_[0] = tomega().name();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find an omega or epsilon field." << nl
            << "buoyancyTurbSource needs omega- or epsilon-based model."
            << exit(FatalError);
    }

    fieldNames_[1] = turbPtr->k()().name();
    applied_.setSize(fieldNames_.size(), false);

    Log << "    Applying buoyancyTurbSource to: "
        << fieldNames_[0] << " and " << fieldNames_[1]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::buoyancyTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(eqn);
        return;
    }

    calcB();

    buoyancyTurbSourceOmega(eqn);
}


void Foam::fv::buoyancyTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(geometricOneField(), rho, eqn, fieldi);
        return;
    }
}


void Foam::fv::buoyancyTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        buoyancyTurbSourceK(alpha, rho, eqn, fieldi);
        return;
    }
}


// ************************************************************************* //
