/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "treeTurbulence.H"
#include "volFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::treeTurbulence::kSource
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvScalarMatrix& eqn,
    const label fieldi
) const
{
    const auto* turbPtr = mesh_.findObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    const volScalarField& k = turbPtr->k();
    const volVectorField& U = turbPtr->U();
    const volScalarField magU(mag(U));
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    // (BSG:Eq. 8)
    eqn +=
        alpha*rho*LAD*Cd*betaP_*pow3(magU)
      - fvm::Sp(alpha*rho*LAD*Cd*betaD_*magU, k);
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::treeTurbulence::epsilonSource
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvScalarMatrix& eqn,
    const label fieldi
) const
{
    const auto* turbPtr = mesh_.findObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    const volScalarField& epsilon = turbPtr->epsilon();
    const volScalarField& k = turbPtr->k();
    const volVectorField& U = turbPtr->U();
    const volScalarField magU(mag(U));
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    // (BSG:Eq. 9)
    eqn +=
        fvm::Sp
        (
            alpha*rho*LAD*Cd*(Ceps1_*betaP_/k*pow3(magU) - Ceps2_*betaD_*magU),
            epsilon
        );
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::treeTurbulence::omegaSource
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvScalarMatrix& eqn,
    const label fieldi
) const
{
    const auto* turbPtr = mesh_.findObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    const volScalarField& omega = turbPtr->omega();
    const volScalarField& k = turbPtr->k();
    const volVectorField& U = turbPtr->U();
    const volScalarField magU(mag(U));
    const volScalarField& Cd = getOrReadField(CdName_);
    const volScalarField& LAD = getOrReadField(LADname_);

    // (derived from BSG:Eq. 9 by assuming epsilon = betaStar*omega*k)
    eqn +=
        fvm::Sp
        (
            alpha*rho*LAD*Cd*betaStar_
           *(Ceps1_*betaP_*pow3(magU) - Ceps2_*betaD_*k*magU),
            omega
        );
}


// ************************************************************************* //
