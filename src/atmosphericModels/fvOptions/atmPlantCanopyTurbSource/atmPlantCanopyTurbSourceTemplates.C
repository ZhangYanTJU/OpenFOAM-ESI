/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "atmPlantCanopyTurbSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::atmPlantCanopyTurbSource::atmPlantCanopyTurbSourceEpsilon
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tepsilon(turbPtr->epsilon());
    const volScalarField& epsilon = tepsilon();
    const volVectorField::Internal& U = turbPtr->U()();

    eqn -= fvm::Sp(alpha()*rho()*(C1_ - C2_)*calcPlantCanopyTerm(U), epsilon);
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::atmPlantCanopyTurbSource::atmPlantCanopyTurbSourceOmega
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
) const
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
    const tmp<volScalarField> tomega(turbPtr->omega());
    const volScalarField& omega = tomega();
    const volVectorField::Internal& U = turbPtr->U()();
    const volScalarField::Internal& gamma =
        mesh_.lookupObjectRef<volScalarField::Internal>
        (
            IOobject::scopedName(turbPtr->type(), "gamma")
        );
    const volScalarField::Internal& beta =
        mesh_.lookupObjectRef<volScalarField::Internal>
        (
            IOobject::scopedName(turbPtr->type(), "beta")
        );

    eqn -= fvm::Sp(alpha()*rho()*(gamma - beta)*calcPlantCanopyTerm(U), omega);
}


// ************************************************************************* //
