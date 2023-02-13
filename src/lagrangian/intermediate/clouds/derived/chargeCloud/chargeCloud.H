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

Class
    Foam::chargeCloud

Description
    Cloud class to introduce multi-phase reacting charged parcels

\*---------------------------------------------------------------------------*/

#ifndef Foam_chargeCloud_H
#define Foam_chargeCloud_H

#include "BaseCloud.H"
#include "KinematicCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "ReactingMultiphaseCloud.H"
#include "ChargeCloud.H"
#include "ParcelCloudModel.H"
#include "chargeParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ParcelCloudModel
    <
        ChargeCloud
        <
            ReactingMultiphaseCloud
            <
                ReactingCloud
                <
                    ThermoCloud
                    <
                        KinematicCloud
                        <
                            BaseCloud
                            <
                                chargeParcel
                            >
                        >
                    >
                >
            >
        >
    > chargeCloud;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
