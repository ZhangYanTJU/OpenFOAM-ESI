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

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::Pstream::broadcast
(
    Type& value,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (is_contiguous_v<Type>)
    {
        UPstream::broadcast
        (
            reinterpret_cast<char*>(&value),
            sizeof(Type),
            communicator
        );
    }
    else
    {
        if (UPstream::master(communicator))
        {
            OPBstream::send(value, communicator);
        }
        else
        {
            IPBstream::recv(value, communicator);
        }
    }
}


template<class Type, unsigned N>
void Foam::Pstream::broadcast
(
    FixedList<Type, N>& list,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (is_contiguous_v<Type>)
    {
        // Size is known and identical on all ranks
        UPstream::broadcast(list.data(), list.size(), communicator);
    }
    else
    {
        // Non-contiguous content - serialize it
        if (UPstream::master(communicator))
        {
            OPBstream::send(list, communicator);
        }
        else
        {
            IPBstream::recv(list, communicator);
        }
    }
}


template<class Type, class... Args>
void Foam::Pstream::broadcasts
(
    const int communicator,
    Type& value,
    Args&&... values
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (!sizeof...(values) && is_contiguous_v<Type>)
    {
        // A single-value and contiguous
        UPstream::broadcast(&value, 1, communicator);
    }
    else
    {
        // Non-contiguous data, or multiple data - needs serialization

        if (UPstream::master(communicator))
        {
            OPBstream::sends
            (
                communicator,
                value,
                std::forward<Args>(values)...
            );
        }
        else
        {
            IPBstream::recvs
            (
                communicator,
                value,
                std::forward<Args>(values)...
            );
        }
    }
}


template<class ListType>
void Foam::Pstream::broadcastList
(
    ListType& list,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        return;
    }
    else if constexpr (is_contiguous_v<typename ListType::value_type>)
    {
        // List data are contiguous
        // 1. broadcast the size
        // 2. resize for receiver list
        // 3. broadcast contiguous contents

        label len(list.size());

        UPstream::mpi_broadcast
        (
            reinterpret_cast<char*>(&len),
            sizeof(label),
            UPstream::dataTypes::type_byte,
            communicator
        );

        if (len)
        {
            // Only broadcast non-empty content
            UPstream::broadcast(list.data(), list.size(), communicator);
        }
    }
    else
    {
        // List data are non-contiguous - serialize/de-serialize

        if (UPstream::master(communicator))
        {
            if (list.empty())
            {
                // Do not serialize if empty.
                // Just broadcast zero-size in a form that IPBstream can expect
                OPBstream::send(Foam::zero{}, communicator);
            }
            else
            {
                OPBstream::send(list, communicator);
            }
        }
        else
        {
            IPBstream is(communicator);
            if (is.remaining() > 0)  // Received a non-empty buffer
            {
                is >> list;
            }
            else
            {
                list.clear();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience wrappers - defined after all specialisations are known

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a broadcasted value (uses a copy internally)
template<class Type>
Type returnBroadcast
(
    const Type& value,
    const int communicator = UPstream::worldComm
)
{
    Type work(value);
    Pstream::broadcast(work, communicator);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
