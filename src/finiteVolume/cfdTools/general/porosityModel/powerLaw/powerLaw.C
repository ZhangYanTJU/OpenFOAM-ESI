/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "powerLaw.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(powerLaw, 0);
        addToRunTimeSelectionTable(porosityModel, powerLaw, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::powerLaw::powerLaw
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    C0_(coeffs_.get<scalar>("C0")),
    C1_(coeffs_.get<scalar>("C1")),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::powerLaw::calcTransformModelData()
{
    // nothing to be transformed
}


void Foam::porosityModels::powerLaw::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(V, rho, U, Udiag);

    force += Udiag*U;
}


void Foam::porosityModels::powerLaw::correct
(
    const volVectorField& U,
    scalarField& Udiag,
    vectorField& Usource,
    bool compressible
) const
{
    const scalarField& V = mesh_.V();

    if (compressible)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(V, rho, U, Udiag);
    }
    else
    {
        apply(V, geometricOneField(), U, Udiag);
    }
}


void Foam::porosityModels::powerLaw::correct
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    scalarField& Udiag,
    vectorField& Usource
) const
{
    const scalarField& V = mesh_.V();

    apply(V, rho, U, Udiag);
}


void Foam::porosityModels::powerLaw::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(rhoName_, U.group())
        );

        apply(rho, U, AU);
    }
    else
    {
        apply(geometricOneField(), U, AU);
    }
}


bool Foam::porosityModels::powerLaw::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);

    return true;
}


// ************************************************************************* //
