/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "pointConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    const char* const pTraits<pointConstraint>::typeName = "pointConstraint";

    defineCompoundTypeName(List<pointConstraint>, pointConstraintList);
    addCompoundToRunTimeSelectionTable
    (
        List<pointConstraint>,
        pointConstraintList
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Binary reading of <label, vector> tuple from a different precision
//
// NOTE: cannot simply use readRawLabel, readRawScalar since the tuple
// may include trailing padding and/or padding between its members

template<class IntType, class FloatType>
static void reading
(
    Foam::Istream& is,
    Foam::pointConstraint* data,
    size_t nElem
)
{
    using namespace Foam;

    // (label, vector) in the specified read precision
    typedef Tuple2<IntType, Vector<FloatType>> inputType;

    for (const pointConstraint* last = data + nElem; data != last; ++data)
    {
        inputType tup;

        is.readRaw(reinterpret_cast<char*>(&tup), sizeof(inputType));

        if constexpr (sizeof(Foam::label) < sizeof(IntType))
        {
            // Narrowing: currently only need (int32_t <- int64_t)
            data->first() = Foam::narrowInt32(tup.first());
        }
        else
        {
            data->first() = tup.first();
        }

        if constexpr (sizeof(Foam::scalar) < sizeof(FloatType))
        {
            // Narrowing: currently only need (float <- double)
            data->second().x() = Foam::narrowFloat(tup.second().x());
            data->second().y() = Foam::narrowFloat(tup.second().y());
            data->second().z() = Foam::narrowFloat(tup.second().z());
        }
        else
        {
            data->second().x() = tup.second().x();
            data->second().y() = tup.second().y();
            data->second().z() = tup.second().z();
        }
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

// Binary reading of Tuple2<label, vector>

template<>
void Foam::Detail::readContiguous<Foam::pointConstraint>
(
    Istream& is,
    char* byteData,
    std::streamsize byteCount
)
{
    if (is.checkNativeSizes())
    {
        // Native label/scalar sizes
        is.read(byteData, byteCount);
    }
    else
    {
        // Non-native label/scalar size
        is.beginRawRead();

        auto* data = reinterpret_cast<pointConstraint*>(byteData);
        const auto nElem = (byteCount/sizeof(pointConstraint));

        switch ((is.labelByteSize() << 8) | is.scalarByteSize())
        {
            case ((sizeof(int32_t) << 8) | sizeof(float)) :
            {
                reading<int32_t, float>(is, data, nElem);
                break;
            }
            case ((sizeof(int32_t) << 8) | sizeof(double)) :
            {
                reading<int32_t, double>(is, data, nElem);
                break;
            }
            case ((sizeof(int64_t) << 8) | sizeof(float)) :
            {
                reading<int32_t, float>(is, data, nElem);
                break;
            }
            case ((sizeof(int64_t) << 8) | sizeof(double)) :
            {
                reading<int32_t, double>(is, data, nElem);
                break;
            }
            default:
            {
                // Unknown combination
                is.fatalCheckNativeSizes(FUNCTION_NAME);
                break;
            }
        }

        is.endRawRead();
    }
}


// ************************************************************************* //
