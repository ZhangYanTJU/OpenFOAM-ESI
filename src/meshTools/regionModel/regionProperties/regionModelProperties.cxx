/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "regionProperties.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// The enable/disable flag for allowFaModels()
static bool allowFaModels_ = true;

bool Foam::regionModels::allowFaModels() noexcept
{
    return allowFaModels_;
}


bool Foam::regionModels::allowFaModels(bool on) noexcept
{
    bool old(allowFaModels_);
    allowFaModels_ = on;
    return old;
}


// ************************************************************************* //
