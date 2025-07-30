/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "HashTableCore.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HashTableCore, 0);
}


// file-scope:
// Minimum internal table size (must be a power of two!)
constexpr int32_t minTableSize = 8;


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::HashTableCore::canonicalSize(const label size) noexcept
{
    // Enforce power of two for fast modulus in hash index calculations.
    // Use unsigned for these calculations.
    //
    // - The lower limit (8) is somewhat arbitrary, but if the hash table
    //   is too small, there will be many direct table collisions.
    // - The upper limit (approx. INT32_MAX/4) must be a power of two,
    //   need not be extremely large for hashing.

    if (size <= minTableSize)
    {
        return (size < 1 ? 0 : minTableSize);
    }
    else if (size > maxTableSize/2)
    {
        return maxTableSize;
    }

    // Determine power-of-two with glibc (may or may not be faster):
    //
    // return (1 << (32-__builtin_clz(int32_t(size-1))));

    if (!(size & (size-1)))
    {
        // Already a power-of-two...
        return size;
    }

    // Non-branching for 32-bit
    // [https://graphics.stanford.edu/~seander/bithacks.html]
    {
        uint32_t n(size);
        --n;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        ++n;
        return n;
    }

    // OLD:
    // Brute-force method
    //
    // uint32_t n(minTableSize);
    // while (n < uint32_t(size))
    // {
    //     n <<= 1;
    // }
    // return n;
}


// ************************************************************************* //
