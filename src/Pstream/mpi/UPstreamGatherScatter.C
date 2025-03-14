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

#include "Pstream.H"
#include "PstreamGlobals.H"
#include "UPstreamWrapping.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static inline bool is_basic_dataType(Foam::UPstream::dataTypes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::dataTypes::Basic_begin)
     && int(id)  < int(Foam::UPstream::dataTypes::Basic_end)
    );
}

namespace
{

using namespace Foam;

// Local function to print some error information
inline void printErrorNonIntrinsic
(
    const char* context,
    UPstream::dataTypes dataTypeId
)
{
    FatalError
        << "Bad input for " << context << ": likely a programming problem\n"
        << "    Non-intrinsic data (" << int(dataTypeId) << ")\n"
        << Foam::endl;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gather
(
    const void* sendData,       // Type checking done by caller
    void* recvData,             // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_gather] :";

        // Appears to be an in-place request
        if
        (
            UPstream::master(communicator)
         && (!sendData || (sendData == recvData))
        )
        {
            Perr<< " (inplace)";
        }

        Perr<< " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular gather

        PstreamDetail::gather
        (
            sendData,
            recvData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_scatter
(
    const void* sendData,       // Type checking done by caller
    void* recvData,             // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_scatter] :";

        // Appears to be an in-place request
        if
        (
            UPstream::master(communicator)
         && (!recvData || (sendData == recvData))
        )
        {
            Perr<< " (inplace)";
        }

        Perr<< " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular scatter

        PstreamDetail::scatter
        (
            sendData,
            recvData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_allgather
(
    void* allData,        // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller

    const int communicator,     // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_allgather] :"
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular all gather

        PstreamDetail::allGather
        (
            allData,
            count,
            datatype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_gatherv
(
    const void* sendData,
    int sendCount,
    void* recvData,
    const UList<int>& recvCounts,
    const UList<int>& recvOffsets,

    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int communicator
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
        )
    )
    {
        FatalErrorInFunction;
        printErrorNonIntrinsic("MPI_Gatherv()", dataTypeId);
        FatalError << Foam::abort(FatalError);
    }

    {
        PstreamDetail::gatherv
        (
            sendData, sendCount,
            recvData, recvCounts, recvOffsets,
            datatype, communicator
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::mpi_scatterv
(
    const void* sendData,
    const UList<int>& sendCounts,
    const UList<int>& sendOffsets,

    void* recvData,
    int recvCount,

    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int communicator
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
        )
    )
    {
        FatalErrorInFunction;
        printErrorNonIntrinsic("MPI_Scatterv()", dataTypeId);
        FatalError << Foam::abort(FatalError);
    }

    {
        PstreamDetail::scatterv
        (
            sendData, sendCounts, sendOffsets,
            recvData, recvCount,
            datatype, communicator
        );
    }
}


// ************************************************************************* //
