/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties()
{}


template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp)
{}


template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict)
{}


template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    eq_(0.0)
{}


template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    eq_(0.0)
{}


template<class ParcelType>
inline Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const scalar eq0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        YGas0,
        YLiquid0,
        YSolid0,
        constProps
    ),
    eq_(eq0)
{}


// * * * * * * * * * * ChargeParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ChargeParcel<ParcelType>::eq() const
{
    return eq_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ChargeParcel<ParcelType>::eq()
{
    return eq_;
}


// ************************************************************************* //
