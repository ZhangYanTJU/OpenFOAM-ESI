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

Class
    Foam::MemoryPool

Description
    Optional memory management using a memory pool such as Umpire
    (https://github.com/LLNL/Umpire).

    When compiled with Umpire, its use can be controlled by the
    \c FOAM_MEMORY_POOL environment variable, or the
    \c memory_pool Optimisation switch (etc/controlDict).

    It currently looks for any of the following entries, in this order:
    - true - same as \em "host"
    - false/none - disabled.
    - \em "host" - uses host memory pool
    - \em "system" - same as \em "host"
    - \em "device" - uses device memory pool
    - \em "managed" - uses managed host/device memory pool
    .

    The parameters "size=nn" and "incr=nn" (in MegaBytes) can be used
    to specify alternatives to the default sizing.

\*---------------------------------------------------------------------------*/

#ifndef Foam_MemoryPool_H
#define Foam_MemoryPool_H

#include <cstdint>  // For size_t

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MemoryPool Declaration
\*---------------------------------------------------------------------------*/

class MemoryPool
{
public:

    // Constructors

        //- Create a memory pool instance (if not already active).
        //  Uses environment or etc/controlDict entry
        static bool create(bool verbose = false);


    // Destructor

        //- Remove the memory pool instance (currently does nothing)
        static void destroy(bool verbose = false);


    // Member Functions

        //- True if pool is active (ie, created and not suspended)
        static bool active() noexcept;

        //- Suspend use of memory pool (for allocation).
        //  \return previous suspend status
        static bool suspend() noexcept;

        //- Resume use of memory pool (if previously active)
        static void resume() noexcept;

        //- Test if given pointer belongs to the pool
        static bool is_pool(void *ptr);

        //- Allocate from pool (if active).
        //  \returns nullptr if the pool is not active
        static void* try_allocate(std::size_t nbytes);

        //- Deallocate a pointer managed by the pool
        //  \returns True if a nullptr (no-op) or when the pointer was
        //      managed by the pool.
        static bool try_deallocate(void *ptr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
