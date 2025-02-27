/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#include "int16.H"
#include "int32.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<int16_t>::typeName = "int16";
const char* const Foam::pTraits<int16_t>::componentNames[] = { "" };

const int16_t Foam::pTraits<int16_t>::zero = 0;
const int16_t Foam::pTraits<int16_t>::one = 1;
const int16_t Foam::pTraits<int16_t>::min = INT16_MIN;
const int16_t Foam::pTraits<int16_t>::max = INT16_MAX;
const int16_t Foam::pTraits<int16_t>::rootMin = INT16_MIN;
const int16_t Foam::pTraits<int16_t>::rootMax = INT16_MAX;


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, int16_t& val)
{
    int32_t parsed;
    is >> parsed;

    val = int16_t(parsed); // narrow
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int16_t val)
{
    return (os << int32_t(val)); // widen
}


// ************************************************************************* //
