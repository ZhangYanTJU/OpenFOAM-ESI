/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "volFields.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    const scalarField& V,
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField& U,
    scalarField& Udiag
) const
{
    const auto& T = mesh_.lookupObject<volScalarField>
    (
        IOobject::groupName(TName_, U.group())
    );

    for (const label zonei : cellZoneIDs_)
    {
        const labelList& cells = mesh_.cellZones()[zonei];

        for (const label celli : cells)
        {
            Udiag[celli] +=
                V[celli]*alpha[celli]*rho[celli]*D_->value(T[celli]);
        }
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField& U,
    tensorField& AU
) const
{
    const auto& T = mesh_.lookupObject<volScalarField>
    (
        IOobject::groupName(TName_, U.group())
    );

    for (const label zonei : cellZoneIDs_)
    {
        const labelList& cells = mesh_.cellZones()[zonei];

        for (const label celli : cells)
        {
            AU[celli] +=
                tensor::I*alpha[celli]*rho[celli]*D_->value(T[celli]);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    const scalarField& V,
    const RhoFieldType& rho,
    const volVectorField& U,
    scalarField& Udiag
) const
{
    if (alphaName_ == "none")
    {
        return apply(V, geometricOneField(), rho, U, Udiag);
    }
    else
    {
        const auto& alpha = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(alphaName_, U.group())
        );

        return apply(V, alpha, rho, U, Udiag);
    }
}


template<class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    const RhoFieldType& rho,
    const volVectorField& U,
    tensorField& AU
) const
{
    if (alphaName_ == "none")
    {
        return apply(geometricOneField(), rho, U, AU);
    }
    else
    {
        const auto& alpha = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(alphaName_, U.group())
        );

        return apply(alpha, rho, U, AU);
    }
}


// ************************************************************************* //
