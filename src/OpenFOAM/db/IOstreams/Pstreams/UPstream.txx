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


template<class Type>
void Foam::UPstream::mpiGather
(
    const Type* sendData,
    Type* recvData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_rank(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else if (!UPstream::is_parallel(communicator))
    {
        // Perform any fallback copying here, while we still know the Type
        if (sendData && recvData && (sendData != recvData))
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_gather
        (
            sendData,   // The data or cmpt pointer
            recvData,   // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class Type>
void Foam::UPstream::mpiScatter
(
    const Type* sendData,
    Type* recvData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_rank(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else if (!UPstream::is_parallel(communicator))
    {
        // Perform any fallback copying here, while we still know the Type
        if (sendData && recvData && (sendData != recvData))
        {
            std::memmove(recvData, sendData, count*sizeof(Type));
        }
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_scatter
        (
            sendData,   // The data or cmpt pointer
            recvData,   // The data or cmpt pointer
            UPstream_dataType<Type>::size(count),
            UPstream_dataType<Type>::datatype_id,
            communicator
        );
    }
}


template<class Type>
void Foam::UPstream::mpiAllGather
(
    Type* allData,
    int count,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing sensible to do
        return;
    }
    else if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        UPstream::mpi_allgather
        (
            allData,    // The data or cmpt pointer
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
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(comm) as well?
        List<T> allValues(1);
        allValues[0] = localValue;
        return allValues;
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot all-gather values for non-contiguous types"
            << " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
        return List<T>();
    }
    else
    {
        // Standard gather with direct MPI communication
        List<T> allValues;

        allValues.resize(UPstream::nProcs(communicator));
        allValues[UPstream::myProcNo(communicator)] = localValue;

        UPstream::mpiAllGather(allValues.data(), 1, communicator);
        return allValues;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T> Foam::UPstream::listGatherValues
(
    const T& localValue,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?
        List<T> allValues(1);
        allValues[0] = localValue;
        return allValues;
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot gather values for non-contiguous types"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);
        return List<T>();
    }
    else
    {
        // Local sizes are identical, can use MPI_Gather
        List<T> allValues;

        if (UPstream::master(communicator))
        {
            allValues.resize(UPstream::nProcs(communicator));
        }

        UPstream::mpiGather(&localValue, allValues.data(), 1, communicator);
        return allValues;
    }
}


template<class T>
T Foam::UPstream::listScatterValues
(
    const UList<T>& allValues,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return own value
        // TBD: only when UPstream::is_rank(communicator) as well?

        if (!allValues.empty())
        {
            return allValues[0];
        }

        return T{};  // Fallback value
    }
    else if constexpr (!is_contiguous_v<T>)
    {
        FatalErrorInFunction
            << "Cannot scatter non-contiguous values"
               " - consider Pstream variant instead" << endl
            << Foam::abort(FatalError);

        return T{};  // Fallback value
    }
    else
    {
        // Local sizes are identical, can use MPI_Scatter

        const label nProcs = UPstream::nProcs(communicator);

        if
        (
            FOAM_UNLIKELY
            (
                UPstream::master(communicator)
             && allValues.size() < nProcs
            )
        )
        {
            FatalErrorInFunction
                << "Attempting to send " << allValues.size()
                << " values to " << nProcs << " processors" << endl
                << Foam::abort(FatalError);
        }

        T localValue{};
        UPstream::mpiScatter(allValues.cdata(), &localValue, 1, communicator);
        return localValue;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::UPstream::mpiGatherv
(
    const Type* sendData,
    int sendCount,
    Type* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        if constexpr (is_contiguous_v<Type>)
        {
            if (sendData && recvData && (sendData != recvData))
            {
                // recvCounts[0] may be invalid - use sendCount instead
                std::memmove(recvData, sendData, sendCount*sizeof(Type));
            }
        }
        // Nothing further to do
    }
    else if constexpr (UPstream_basic_dataType<Type>::value)
    {
        // Restrict to basic (or aliased) MPI types to avoid recalculating
        // the list of counts/offsets.

        UPstream::mpi_gatherv
        (
            sendData,
            sendCount,
            recvData,
            recvCounts,
            recvOffsets,

            UPstream_basic_dataType<Type>::datatype_id,
            communicator
        );
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<Type>, "Only basic MPI data types"
        );
    }
}


template<class Type>
void Foam::UPstream::mpiScatterv
(
    const Type* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,
    Type* recvData,
    int recvCount,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        if constexpr (is_contiguous_v<Type>)
        {
            if (sendData && recvData && (sendData != recvData))
            {
                std::memmove(recvData, sendData, recvCount*sizeof(Type));
            }
        }
        // Nothing further to do
    }
    else if constexpr (UPstream_basic_dataType<Type>::value)
    {
        // Restrict to basic (or aliased) MPI types to avoid recalculating
        // the list of counts/offsets.

        UPstream::mpi_scatterv
        (
            sendData,
            sendCounts,
            sendOffsets,
            recvData,
            recvCount,

            UPstream_basic_dataType<Type>::datatype_id,
            communicator
        );
    }
    else
    {
        static_assert
        (
            stdFoam::dependent_false_v<Type>, "Only basic MPI data types"
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void Foam::UPstream::mpiReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_reduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_allreduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator
        );
    }
}


template<class T>
void Foam::UPstream::mpiAllReduce
(
    T values[],
    int count,
    const UPstream::opCodes opCodeId,
    const int communicator,
    UPstream::Request& req
)
{
    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // Use element or component type (or byte-wise) for data type
        // Restricted to basic data types
        UPstream::mpi_allreduce
        (
            values,     // The data or cmpt pointer
            UPstream_basic_dataType<T>::size(count),
            UPstream_basic_dataType<T>::datatype_id,
            opCodeId,
            communicator,
           &req
        );
    }
}


// ************************************************************************* //
