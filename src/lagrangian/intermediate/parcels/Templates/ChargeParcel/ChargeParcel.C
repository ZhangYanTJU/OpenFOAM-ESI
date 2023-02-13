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

#include "ChargeParcel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ChargeParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ChargeParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ChargeParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    ParcelType::calc(cloud, td, dt);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const ChargeParcel<ParcelType>& p
)
:
    ParcelType(p),
    eq_(p.eq_)
{}


template<class ParcelType>
Foam::ChargeParcel<ParcelType>::ChargeParcel
(
    const ChargeParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    eq_(p.eq_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ChargeParcelIO.C"

// ************************************************************************* //
