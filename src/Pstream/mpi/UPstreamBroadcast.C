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

#include "UPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

bool Foam::UPstream::mpi_broadcast
(
    void* buf,                      // Type checking done by caller
    std::streamsize count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const int communicator          // Index into MPICommunicators_
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return true;
    }

    //Needed?  PstreamGlobals::checkCommunicator(communicator, 0);

    // Without MPI_Bcast_c()
    if (FOAM_UNLIKELY(count > std::streamsize(INT_MAX)))
    {
        FatalErrorInFunction
            << "Broadcast size " << label(count)
            << " exceeds INT_MAX bytes" << Foam::endl
            << Foam::abort(FatalError);
        return false;
    }

    const bool withTopo =
    (
        UPstream::usingTopoControl(UPstream::topoControls::broadcast)
     && UPstream::usingNodeComms(communicator)
    );

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_broadcast] :"
            << " type:" << int(dataTypeId)
            << " count:" << label(count)
            << " comm:" << communicator
            << " topo:" << withTopo << Foam::endl;
    }

    int returnCode = MPI_SUCCESS;

    profilingPstream::beginTiming();

    if (withTopo)
    {
        // Topological broadcast

        for
        (
            const int subComm :
            // std::initializer_list<int>
            {
                UPstream::commInterNode_,   // Stage 1: between nodes
                UPstream::commLocalNode_    // Stage 2: within a node
            }
        )
        {
            if (UPstream::is_parallel(subComm))
            {
                if (FOAM_UNLIKELY(UPstream::debug))
                {
                    Perr<< "[mpi_broadcast] :"
                        << " type:" << int(dataTypeId)
                        << " count:" << label(count)
                        << " comm:" << subComm
                        << " substage" << Foam::endl;
                }

                returnCode = MPI_Bcast
                (
                    buf,
                    count,
                    datatype,
                    0,  // (root rank) == UPstream::masterNo()
                    PstreamGlobals::MPICommunicators_[subComm]
                );
            }
        }
    }
    else
    {
        // Regular broadcast
        // OR: PstreamDetail::broadcast(buf, count, datatype, communicator);

        returnCode = MPI_Bcast
        (
            buf,
            count,
            datatype,
            0,  // (root rank) == UPstream::masterNo()
            PstreamGlobals::MPICommunicators_[communicator]
        );
    }

    profilingPstream::addBroadcastTime();

    return (returnCode == MPI_SUCCESS);
}


// ************************************************************************* //
