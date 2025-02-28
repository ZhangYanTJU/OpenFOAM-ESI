/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::UPstream::broadcast
(
    Type* buffer,
    std::streamsize count,
    const int communicator
)
{
    // Likely no reason to check for nullptr
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return true;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        // Also report parameters to silence compiler warnings about unused
        FatalErrorInFunction
            << "Invalid for non-contiguous data types."
            << " buffer:" << (buffer != nullptr)
            << " count:" << count
            << Foam::abort(FatalError);
        return false;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        return UPstream::mpi_broadcast
        (
            buffer,     // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class T>
Foam::List<T> Foam::UPstream::allGatherValues
(
    const T& localValue,
    const label comm
)
{
    if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot all-gather values for non-contiguous types"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
    }

    List<T> allValues;

    if (UPstream::is_parallel(comm))
    {
        allValues.resize(UPstream::nProcs(comm));
        allValues[UPstream::myProcNo(comm)] = localValue;

        UPstream::mpiAllGather(allValues.data_bytes(), sizeof(T), comm);
    }
    else
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(comm) as well?
        allValues.resize(1);
        allValues[0] = localValue;
    }

    return allValues;
}


template<class T>
Foam::List<T> Foam::UPstream::listGatherValues
(
    const T& localValue,
    const label comm
)
{
    if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot gather values for non-contiguous types"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
    }

    List<T> allValues;

    if (UPstream::is_parallel(comm))
    {
        if (UPstream::master(comm))
        {
            allValues.resize(UPstream::nProcs(comm));
        }

        UPstream::mpiGather
        (
            reinterpret_cast<const char*>(&localValue),
            allValues.data_bytes(),
            sizeof(T),  // The send/recv size per rank
            comm
        );
    }
    else
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(comm) as well?
        allValues.resize(1);
        allValues[0] = localValue;
    }

    return allValues;
}


template<class T>
T Foam::UPstream::listScatterValues
(
    const UList<T>& allValues,
    const label comm
)
{
    if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot scatter values for non-contiguous types"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
    }

    T localValue{};

    if (UPstream::is_parallel(comm))
    {
        const label numProc = UPstream::nProcs(comm);

        if (UPstream::master(comm) && allValues.size() < numProc)
        {
            FatalErrorInFunction
                << "Attempting to send " << allValues.size()
                << " values to " << numProc << " processors" << endl
                << Foam::abort(FatalError);
        }

        UPstream::mpiScatter
        (
            allValues.cdata_bytes(),
            reinterpret_cast<char*>(&localValue),
            sizeof(T),  // The send/recv size per rank
            comm
        );
    }
    else
    {
        // non-parallel: return first value
        // TBD: only when UPstream::is_rank(comm) as well?

        if (!allValues.empty())
        {
            return allValues[0];
        }
     }

     return localValue;
}


// ************************************************************************* //
