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

Application
    Test-globalIndex3

Description
    Tests for globalIndex with node-wise splitting

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "globalMeshData.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "Random.H"
#include "openfoam_mpi.H"

// pre-scan for "-split-size NUM"
int option_splitsize(int argc, char *argv[])
{
    int ivalue = -1;
    for (int argi = 1; argi < argc-1; ++argi)
    {
        if (strcmp(argv[argi], "-split-size") == 0)
        {
            ++argi;
            ivalue = atoi(argv[argi]);
        }
    }

    return ivalue;
}

using namespace Foam;

template<class T>
void printList(Ostream& os, const UList<T>& list)
{
    os << list.size() << " " << flatOutput(list) << nl;
}

void printGlobalIndex(Ostream& os, const globalIndex& gi)
{
    printList(os, gi.offsets());
}


bool UPstream_usingHostComms(const label communicator = UPstream::worldComm)
{
    // Starting point must be "real" world-communicator
    // ("real" means without any local trickery with worldComm)
    // Avoid corner cases:
    // - everthing is on one node
    // - everthing is on different nodes

    label numNodes_ = UPstream::nProcs(UPstream::commInterHost());
    UPstream::broadcast
    (
        reinterpret_cast<char*>(&numNodes_),
        sizeof(label),
        UPstream::commWorld()
    );

    return
    (
        // (UPstream::hostCommsEnabled_ > 0)
        UPstream::parRun()
     && (UPstream::commWorld() == communicator)
        // More than one node and some processes do share nodes
     && (numNodes_ > 1)  // ...
    );
}


template<class ProcIDsContainer, class Type>
void globalIndexGather
(
    const labelUList& off,  // needed on master only
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag,
    const UPstream::commsTypes preferredCommsType,
    bool useWindow = false
)
{
    // low-level: no parRun guard

    // Protection for disjoint calls
    if (FOAM_UNLIKELY(!UPstream::is_rank(comm)))
    {
        FatalErrorInFunction
            << "Calling with process not on the communicator"
            << Foam::abort(FatalError);
    }

    // Cannot use non-blocking for non-contiguous data.
    const UPstream::commsTypes commsType =
    (
        (
            !is_contiguous_v<Type>
         && UPstream::commsTypes::nonBlocking == preferredCommsType
        )
      ? UPstream::commsTypes::scheduled
      : preferredCommsType
    );


    const label startOfRequests = UPstream::nRequests();

    const int masterProci = procIDs.size() ? procIDs[0] : 0;


    // Very hard-coded at the moment
    int returnCode = MPI_SUCCESS;
    const int nCmpts = pTraits<Type>::nComponents;

    MPI_Win win;
    MPI_Datatype dataType = MPI_DOUBLE;
    if (useWindow)
    {
        using cmptType = typename pTraits<Type>::cmptType;

        if (std::is_same<float, cmptType>::value)
        {
            dataType = MPI_FLOAT;
        }
        else if (std::is_same<double, cmptType>::value)
        {
            dataType = MPI_DOUBLE;
        }
        else
        {
            // Not supported
            useWindow = false;
        }
    }

    if (useWindow)
    {
        MPI_Comm mpiComm =
            PstreamUtils::Cast::to_mpi(UPstream::Communicator::lookup(comm));

        char commName[MPI_MAX_OBJECT_NAME];
        int nameLen = 0;

        if
        (
            MPI_COMM_NULL != mpiComm
         && MPI_SUCCESS == MPI_Comm_get_name(mpiComm, commName, &nameLen)
         && (nameLen > 0)
        )
        {
            Pout<< "window on " << commName << nl;
        }

        if (UPstream::myProcNo(comm) == masterProci || fld.empty())
        {
            // Collective
            returnCode = MPI_Win_create
            (
                nullptr,
                0,
                1,  // disp_units
                MPI_INFO_NULL,
                mpiComm,
               &win
            );
        }
        else
        {
            // Collective
            returnCode = MPI_Win_create
            (
                const_cast<char *>(fld.cdata_bytes()),
                fld.size_bytes(),
                sizeof(Type),  // disp_units
                MPI_INFO_NULL,
                mpiComm,
               &win
            );
        }

        if (MPI_SUCCESS != returnCode || MPI_WIN_NULL == win)
        {
            FatalErrorInFunction
                << "MPI_Win_create() failed"
                << Foam::abort(FatalError);
            // return nullptr;
        }
    }


    if (UPstream::myProcNo(comm) == masterProci)
    {
        allFld.resize_nocopy(off.back());  // == totalSize()

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.

        SubList<Type>(allFld, off[1]-off[0], off[0]) =
            SubList<Type>(fld, off[1]-off[0]);

        if (useWindow)
        {
            MPI_Win_lock_all(MPI_MODE_NOCHECK, win);
        }

        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

            if (procSlot.empty())
            {
                // Nothing to do
            }
            else if (useWindow)
            {
                returnCode = MPI_Get
                (
                    // origin
                    procSlot.data(),
                    procSlot.size()*(nCmpts),
                    dataType,

                    // target
                    i,
                    0,
                    procSlot.size()*(nCmpts),
                    dataType,
                    win
                );

                if (MPI_SUCCESS != returnCode)
                {
                    FatalErrorInFunction
                        << "MPI_Get failed"
                        << Foam::abort(FatalError);
                    // return nullptr;
                }
            }
            else if constexpr (is_contiguous_v<Type>)
            {
                UIPstream::read
                (
                    commsType,
                    procIDs[i],
                    procSlot,
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(procSlot, procIDs[i], tag, comm);
                IPstream::recv(procSlot, i, tag, comm);
            }
        }

        if (useWindow)
        {
            MPI_Win_unlock_all(win);
        }
    }
    else if (!useWindow)
    {
        if (fld.empty())
        {
            // Nothing to do
        }
        else if constexpr (is_contiguous_v<Type>)
        {
            UOPstream::write
            (
                commsType,
                masterProci,
                fld,
                tag,
                comm
            );
        }
        else
        {
            OPstream::send(fld, commsType, masterProci, tag, comm);
        }
    }

    if (useWindow)
    {
        // Collective
        MPI_Win_free(&win);
    }

    if (commsType == UPstream::commsTypes::nonBlocking)
    {
        // Wait for outstanding requests
        UPstream::waitRequests(startOfRequests);
    }
}


// Populate the inter-node/intra-node offsets
static void splitOffsets
(
    const globalIndex& gi,

    // Offsets between node leaders
    labelList& interOffsets,

    // Offsets with a node
    labelList& intraOffsets
)
{
    const auto interComm = UPstream::commInterHost();

    // Only spread information to the node leaders
    if (UPstream::is_rank(interComm))
    {
        const auto& procIds = UPstream::procID(interComm);
        const auto myInterRank = UPstream::myProcNo(interComm);
        const auto numInter = UPstream::nProcs(interComm);

        labelList allOffsets;
        if (UPstream::master(interComm))
        {
            allOffsets = offsets;
        }

        // FUTURE: will have UPstream::numNodes() everywhere
        // and can size allOffsets accordingly - avoids double broadcast

        Pstream::broadcastList(allOffsets, interComm);


        // intra-nodes offsets:
        if (!allOffsets.empty())
        {
            if (myInterRank+1 < numInter)
            {
                intraOffsets = allOffsets.slice
                (
                    procIds[myInterRank],
                    (procIds[myInterRank+1] - procIds[myInterRank] + 1)
                );
            }
            else
            {
                intraOffsets = allOffsets.slice(procIds[myInterRank]);
            }
        }
        if (!intraOffsets.empty())
        {
            const auto start0 = intraOffsets.front();
            for (auto& val : intraOffsets)
            {
                val -= start0;
            }
        }

        // inter-nodes offsets:
        {
            interOffsets.resize_nocopy(procIds.size()+1);

            forAll(procIds, i)
            {
                interOffsets[i] = allOffsets[procIds[i]];
            }
            interOffsets.back() = allOffsets.back();
        }
    }
    else
    {
        interOffsets.clear();
        intraOffsets.clear();
    }
}


template<class Type>
void globalIndexGather
(
    const globalIndex& gi,
    const UList<Type>& sendData,
    List<Type>& allData,
    const int tag,
    const UPstream::commsTypes commsType,
    const label comm = UPstream::worldComm,
    bool useWindow = false
)
{
    if (!UPstream::parRun())
    {
        // Serial: direct copy
        allData = sendData;
        return;
    }

    const auto& offsets = gi.offsets();  // needed on master only

    const bool usingHostComms = UPstream_usingHostComms(comm);

    Info<< "Using host-comms: " << usingHostComms << nl;

    if (UPstream_usingHostComms(comm))
    {
        // Stage 0 : The inter-node/intra-node offsets
        labelList interOffsets;
        labelList intraOffsets;

        splitOffsets(gi, interOffsets, interOffsets);

        List<Type> nodeData;

        // Stage 1 : Gather intra-node data
        {
            const auto currComm = UPstream::commIntraHost();

            globalIndexGather
            (
                intraOffsets,  // needed on master only
                currComm,
                UPstream::allProcs(currComm),
                sendData,
                nodeData,
                tag,
                commsType,
                useWindow
            );
        }

        // Stage 2 : Gather inter-node data
        if (UPstream::is_rank(UPstream::commInterHost()))
        {
            const auto currComm = UPstream::commInterHost();

            globalIndexGather
            (
                interOffsets,  // needed on master only
                currComm,
                UPstream::allProcs(currComm),
                nodeData,
                allData,
                tag,
                commsType,
                useWindow
            );
        }
    }
    else
    {
        globalIndexGather
        (
            offsets,  // needed on master only
            comm,
            UPstream::allProcs(comm),
            sendData,
            allData,
            tag,
            commsType,
            useWindow
        );
    }

    if (!UPstream::master(comm))
    {
        allData.clear();  // safety: zero-size on non-master
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("Set UPstream::debug level");
    argList::addOption("split-size", "NUM", "split with ncores/node");
    argList::addBoolOption("window", "get data via window");

    // Check -verbose before initialisation
    UPstream::debug = argList::verbose(argc, argv);

    // Check -split-size before initialisation
    {
        int splitSize = option_splitsize(argc, argv);

        if (splitSize >= 0)
        {
            UPstream::hostCommsEnabled_ = splitSize;
        }
    }

    #include "setRootCase.H"

    const bool useLocalComms = UPstream_usingHostComms();
    const bool useWindow = args.found("window");

    Info<< nl
        << "Getting local-comms: " << Switch::name(useLocalComms) << nl
        << "Getting data with window: " << Switch::name(useWindow) << nl
        << nl;


    Random rng(31 + 2*UPstream::myProcNo());

    const label localSize = (5*rng.position<label>(1, 15));

    globalIndex globIndex
    (
        globalIndex::gatherOnly{},
        localSize,
        UPstream::commWorld()
    );

    Info<< "global-index: ";
    printGlobalIndex(Info, globIndex);

    // Splitting offsets
    labelList allOffsets;
    labelList interOffsets;
    labelList intraOffsets;

    labelRange nodeRange;

    if (UPstream::is_rank(UPstream::commInterHost()))
    {
        if (UPstream::master(UPstream::commInterHost()))
        {
            allOffsets = globIndex.offsets();
        }
        Pstream::broadcastList(allOffsets, UPstream::commInterHost());

        const int myInterRank = UPstream::myProcNo(UPstream::commInterHost());
        const auto& procIds = UPstream::procID(UPstream::commInterHost());

        // intra-nodes offsets:
        if (!allOffsets.empty())
        {
            if (myInterRank+1 < UPstream::nProcs(UPstream::commInterHost()))
            {
                intraOffsets = allOffsets.slice
                (
                    procIds[myInterRank],
                    (procIds[myInterRank+1] - procIds[myInterRank] + 1)
                );
            }
            else
            {
                intraOffsets = allOffsets.slice(procIds[myInterRank]);
            }
        }
        if (!intraOffsets.empty())
        {
            const auto start0 = intraOffsets.front();
            for (auto& val : intraOffsets)
            {
                val -= start0;
            }
        }

        // inter-nodes offsets:
        {
            interOffsets.resize(procIds.size()+1);

            forAll(procIds, i)
            {
                interOffsets[i] = allOffsets[procIds[i]];
            }
            interOffsets.back() = allOffsets.back();
        }

        // For reporting
        {
            auto nInterRank = UPstream::nProcs(UPstream::commInterHost());

            nodeRange.reset
            (
                procIds[myInterRank],
                (
                    myInterRank < nInterRank-1
                  ? procIds[myInterRank+1]
                  : UPstream::nProcs(UPstream::commWorld())
                ) - procIds[myInterRank]
            );
        }

    }

    // Pstream::broadcastList(intraOffsets, UPstream::commIntraHost());

    Pout<< "node-range: " << nodeRange << nl;
    Pout<< "all-offset: "; printList(Pout, allOffsets);
    Pout<< "inter-offset: "; printList(Pout, interOffsets);
    Pout<< "intra-offset: "; printList(Pout, intraOffsets);

    Field<scalar> allData;
    Field<scalar> localFld(localSize, scalar(UPstream::myProcNo()));

    // globIndex.gather
    // (
    //     localFld,
    //     allData,
    //     UPstream::msgType(),
    //     UPstream::commsTypes::nonBlocking,
    //     UPstream::commWorld()
    // );

    globalIndexGather
    (
        globIndex,
        localFld,
        allData,
        UPstream::msgType(),
        UPstream::commsTypes::nonBlocking,
        UPstream::commWorld(),
        useWindow
    );

    Pout<< "local: " << flatOutput(localFld) << nl;
    Info<< "field: " << flatOutput(allData) << nl;

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
