/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "uint32.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<uint32_t>::componentNames[] = { "" };

const uint32_t Foam::pTraits<uint32_t>::zero = 0;
const uint32_t Foam::pTraits<uint32_t>::one = 1;
const uint32_t Foam::pTraits<uint32_t>::min = 0;
const uint32_t Foam::pTraits<uint32_t>::max = UINT32_MAX;
const uint32_t Foam::pTraits<uint32_t>::rootMin = 0;
const uint32_t Foam::pTraits<uint32_t>::rootMax = UINT32_MAX;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pTraits<uint32_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
