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
    Foam::UPstream::Window

Description
    An opaque wrapper for MPI_Win with a vendor-independent
    representation and without any \c <mpi.h> header dependency.

Note
    The MPI standard states that MPI_Win is always an opaque object.
    Generally it is either an integer (eg, mpich) or a pointer (eg, openmpi).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamWindow_H
#define Foam_UPstreamWindow_H

#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class UPstream::Window Declaration
\*---------------------------------------------------------------------------*/

class UPstream::Window
{
public:

    // Public Types

        //- Storage for MPI_Win (as integer or pointer)
        typedef std::intptr_t value_type;


private:

    // Private Data

        //- The MPI_Win (as wrapped value)
        value_type value_;

public:

    // Generated Methods

        //- Copy construct
        Window(const Window&) noexcept = default;

        //- Move construct
        Window(Window&&) noexcept = default;

        //- Copy assignment
        Window& operator=(const Window&) noexcept = default;

        //- Move assignment
        Window& operator=(Window&&) noexcept = default;


    // Member Operators

        //- Test for equality
        bool operator==(const Window& rhs) const noexcept
        {
            return (value_ == rhs.value_);
        }

        //- Test for inequality
        bool operator!=(const Window& rhs) const noexcept
        {
            return (value_ != rhs.value_);
        }


    // Constructors

        //- Default construct as MPI_WIN_NULL
        Window() noexcept;

        //- Construct from MPI_Win (as pointer type)
        explicit Window(const void* p) noexcept
        :
            value_(reinterpret_cast<value_type>(p))
        {}

        //- Construct from MPI_Win (as integer type)
        explicit Window(value_type val) noexcept
        :
            value_(val)
        {}


    // Member Functions

    // Basic handling

        //- Return raw value
        value_type value() const noexcept { return value_; }

        //- Return as pointer value
        const void* pointer() const noexcept
        {
            return reinterpret_cast<const void*>(value_);
        }

        //- True if not equal to MPI_WIN_NULL
        bool good() const noexcept;

        //- Reset to default constructed value (MPI_WIN_NULL)
        void reset() noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
