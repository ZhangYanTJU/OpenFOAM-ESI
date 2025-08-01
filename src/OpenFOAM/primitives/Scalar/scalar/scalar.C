/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

#include "scalar.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::scalar Foam::readScalarOrDefault(Istream& is, const scalar defaultValue)
{
    if (is.good())
    {
        token tok(is);

        // NB: does not handle separated '-' (or '+') prefixes
        // like operator>>(Istream&, scalar&) does

        if (tok.isNumber())
        {
            return tok.number();
        }

        is.putBack(tok);
    }

    return defaultValue;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Binary reading with narrowing/widening
template<class readType, class dataType>
static void reading(Foam::Istream& is, dataType* data, size_t nElem)
{
    if constexpr (sizeof(dataType) == sizeof(readType))
    {
        // Read uses the native data size
        is.readRaw(reinterpret_cast<char*>(data), nElem*sizeof(dataType));
    }
    else
    {
        for (const dataType* last = data + nElem; data != last; ++data)
        {
            readType val;
            is.readRaw(reinterpret_cast<char*>(&val), sizeof(readType));

            if constexpr (sizeof(dataType) < sizeof(readType))
            {
                // Narrowing: currently only need (float <- double)
                *data = Foam::narrowFloat(val);
            }
            else
            {
                // Type widening
                *data = dataType(val);
            }
        }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::readRawScalar(Istream& is, scalar* data, size_t nElem)
{
    // No check for binary vs ascii, the caller knows what they are doing

    switch (is.scalarByteSize())
    {
        case sizeof(float):
        {
            reading<float, scalar>(is, data, nElem);
            break;
        }
        case sizeof(double):
        {
            reading<double, scalar>(is, data, nElem);
            break;
        }
        default:
        {
            // Cannot recover from this
            FatalIOErrorInFunction(is)
                << "Currently no code to read float" << (8*is.scalarByteSize())
                << " as float" << (8*sizeof(scalar)) << nl
                << abort(FatalIOError);
        }
    }
}


// ************************************************************************* //
