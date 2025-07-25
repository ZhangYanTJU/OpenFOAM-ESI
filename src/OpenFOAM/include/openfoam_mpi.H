/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

Description
    Header for low-level interfaces between MPI and OpenFOAM.
    The detail interfaces are subject to change.

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamMPI_H
#define Foam_UPstreamMPI_H

#include "UPstream.H"

// Include MPI without any C++ bindings
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamUtils
{

// Casting helpers
struct Cast
{
    //- Cast UPstream::Communicator to MPI_Comm
    template<typename Type = MPI_Comm>
    static Type to_mpi(UPstream::Communicator arg) noexcept
    {
        if constexpr (std::is_pointer_v<Type>)
        {
            return reinterpret_cast<Type>(arg.value());
        }
        else  // std::is_integral_v<Type>
        {
            return static_cast<Type>(arg.value());
        }
    }

    //- Cast UPstream::Request to MPI_Request
    template<typename Type = MPI_Request>
    static Type to_mpi(UPstream::Request arg) noexcept
    {
        if constexpr (std::is_pointer_v<Type>)
        {
            return reinterpret_cast<Type>(arg.value());
        }
        else  // std::is_integral_v<Type>
        {
            return static_cast<Type>(arg.value());
        }
    }

    //- Cast UPstream::Window to MPI_Win
    template<typename Type = MPI_Win>
    static Type to_mpi(UPstream::Window arg) noexcept
    {
        if constexpr (std::is_pointer_v<Type>)
        {
            return reinterpret_cast<Type>(arg.value());
        }
        else  // std::is_integral_v<Type>
        {
            return static_cast<Type>(arg.value());
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
