/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
#include "PstreamReduceOps.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"
#include "int.H"
#include "UPstreamWrapping.H"
#include "collatedFileOperation.H"

#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>

#undef Pstream_use_MPI_Get_count

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The min value and default for MPI buffer length
constexpr int minBufLen = 20000000;

// Track size of attached MPI buffer
static int attachedBufLen = 0;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Attach user-defined send buffer
static void attachOurBuffers()
{
#ifndef SGIMPI
    if (attachedBufLen)
    {
        return;  // Already attached
    }

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str, len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    char* buf = new char[len];

    if (MPI_SUCCESS == MPI_Buffer_attach(buf, len))
    {
        // Properly attached
        attachedBufLen = len;

        if (Foam::UPstream::debug)
        {
            Foam::Perr<< "UPstream::init : buffer-size " << len << '\n';
        }
    }
    else
    {
        delete[] buf;
        Foam::Perr<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


// Remove an existing user-defined send buffer
// IMPORTANT:
//     This operation will block until all messages currently in the
//     buffer have been transmitted.
static void detachOurBuffers()
{
#ifndef SGIMPI
    if (!attachedBufLen)
    {
        return;  // Nothing to detach
    }

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

    char* buf = nullptr;
    int len = 0;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        // This was presumably the buffer that we attached
        // and not someone else.
        delete[] buf;
    }

    // Nothing attached
    attachedBufLen = 0;
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        if (UPstream::debug)
        {
            Perr<< "UPstream::initNull : was already initialized\n";
        }
    }
    else
    {
        // Not already initialized

        MPI_Init_thread
        (
            nullptr,    // argc
            nullptr,    // argv
            MPI_THREAD_SINGLE,
            &flag       // provided_thread_support
        );

        ourMpi = true;
    }

    // Could also attach buffers etc.

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int numprocs = 0, myRank = 0;
    int provided_thread_support = 0;
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized.
        // Warn if we've called twice, but skip if initialized externally

        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already initialized - cannot perform MPI_Init" << nl
                << "This could indicate an application programming error!"
                << endl;

            return true;
        }
        else if (UPstream::debug)
        {
            Perr<< "UPstream::init : was already initialized\n";
        }
    }
    else
    {
        MPI_Init_thread
        (
            &argc,
            &argv,
            (
                needsThread
              ? MPI_THREAD_MULTIPLE
              : MPI_THREAD_SINGLE
            ),
            &provided_thread_support
        );

        ourMpi = true;
    }

    // Check argument list for local world
    label worldIndex = -1;
    word world;
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "-world") == 0)
        {
            worldIndex = argi++;
            if (argi >= argc)
            {
                FatalErrorInFunction
                    << "Missing world name to argument \"world\""
                    << Foam::abort(FatalError);
            }
            world = argv[argi];
            break;
        }
    }

    // Filter 'world' option
    if (worldIndex != -1)
    {
        for (label i = worldIndex+2; i < argc; i++)
        {
            argv[i-2] = argv[i];
        }
        argc -= 2;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (UPstream::debug)
    {
        Perr<< "UPstream::init :"
            << " thread-support : requested:" << needsThread
            << " obtained:"
            << (
                   (provided_thread_support == MPI_THREAD_SINGLE)
                 ? "SINGLE"
                 : (provided_thread_support == MPI_THREAD_SERIALIZED)
                 ? "SERIALIZED"
                 : (provided_thread_support == MPI_THREAD_MULTIPLE)
                 ? "MULTIPLE"
                 : "other"
               )
            << " procs:" << numprocs
            << " rank:" << myRank
            << " world:" << world << endl;
    }

    if (worldIndex == -1 && numprocs <= 1)
    {
        FatalErrorInFunction
            << "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    if (worldIndex != -1)
    {
        // During startup, so commWorld() == commGlobal()

        wordList worlds(numprocs);
        worlds[UPstream::myProcNo(UPstream::commGlobal())] = world;
        Pstream::gatherList
        (
            worlds,
            UPstream::msgType(),
            UPstream::commGlobal()
        );

        // Compact
        if (UPstream::master(UPstream::commGlobal()))
        {
            DynamicList<word> worldNames(numprocs);
            worldIDs_.resize_nocopy(numprocs);

            forAll(worlds, proci)
            {
                const word& world = worlds[proci];

                worldIDs_[proci] = worldNames.find(world);

                if (worldIDs_[proci] == -1)
                {
                    worldIDs_[proci] = worldNames.size();
                    worldNames.push_back(world);
                }
            }

            allWorlds_.transfer(worldNames);
        }
        Pstream::broadcasts(UPstream::commGlobal(), allWorlds_, worldIDs_);

        const label myWorldId =
            worldIDs_[UPstream::myProcNo(UPstream::commGlobal())];

        DynamicList<label> subRanks;
        forAll(worldIDs_, proci)
        {
            if (worldIDs_[proci] == myWorldId)
            {
                subRanks.push_back(proci);
            }
        }

        // Allocate new world-communicator with comm-global as its parent.
        const label subComm =
            UPstream::allocateCommunicator(UPstream::commGlobal(), subRanks);

        // Override world communicator.
        // - the updated (const) world comm does not change after this.
        UPstream::constWorldComm_ = subComm;
        UPstream::worldComm = subComm;

        // For testing: warn use of non-worldComm
        UPstream::warnComm = UPstream::worldComm;

        // MPI_COMM_SELF : the processor number wrt the new world communicator
        if (procIDs_[UPstream::commSelf()].size())
        {
            procIDs_[UPstream::commSelf()].front() =
                UPstream::myProcNo(subComm);
        }

        // Provide some names for these communicators
        if
        (
            MPI_COMM_NULL
         != PstreamGlobals::MPICommunicators_[UPstream::commGlobal()]
        )
        {
            MPI_Comm_set_name
            (
                PstreamGlobals::MPICommunicators_[UPstream::commGlobal()],
                "<global>"
            );
        }
        if
        (
            MPI_COMM_NULL
         != PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_]
        )
        {
            MPI_Comm_set_name
            (
                PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_],
                ("world=" + world).c_str()
            );
        }

        if (UPstream::debug)
        {
            // Check
            const auto mpiNewComm = PstreamGlobals::MPICommunicators_[subComm];
            int subNumProcs, subRank;
            MPI_Comm_size(mpiNewComm, &subNumProcs);
            MPI_Comm_rank(mpiNewComm, &subRank);

            Perr<< "UPstream::init : in world:" << world
                << " using local communicator:" << subComm
                << " rank " << subRank
                << " of " << subNumProcs
                << endl;
        }

        // Override Pout prefix (move to setParRun?)
        Pout.prefix() = '[' + world + '/' +  name(myProcNo(subComm)) + "] ";
        Perr.prefix() = Pout.prefix();
    }
    else
    {
        // All processors use world 0
        worldIDs_.resize_nocopy(numprocs);
        worldIDs_ = 0;

        // Provide some names for these communicators
        if
        (
            MPI_COMM_NULL
         != PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_]
        )
        {
            MPI_Comm_set_name
            (
                PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_],
                "<openfoam:world>"
            );
        }
    }


    // Define inter-host and intra-host communicators.
    if (UPstream::hostCommsEnabled_ >= 8)
    {
        // Debugging: split with given number per node
        setHostCommunicators(UPstream::hostCommsEnabled_);
    }
    #ifndef MSMPI_VER  /* Uncertain if this would work with MSMPI */
    else if (UPstream::hostCommsEnabled_ == 2)
    {
        // Defined based on shared-memory hardware information
        setSharedMemoryCommunicators();
    }
    #endif
    else
    {
        // Defined based on hostname, even if nominally disabled
        setHostCommunicators();
    }


    // Provide some names for these communicators
    if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[interHostComm_])
    {
        MPI_Comm_set_name
        (
            PstreamGlobals::MPICommunicators_[interHostComm_],
            "<openfoam:inter-host>"
        );
    }
    if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[intraHostComm_])
    {
        MPI_Comm_set_name
        (
            PstreamGlobals::MPICommunicators_[intraHostComm_],
            "<openfoam:intra-host>"
        );
    }

    attachOurBuffers();

    return true;
}


void Foam::UPstream::shutdown(int errNo)
{
    int flag = 0;

    MPI_Initialized(&flag);
    if (!flag)
    {
        // MPI not initialized - we have nothing to do
        return;
    }

    MPI_Finalized(&flag);
    if (flag)
    {
        // MPI already finalized - we have nothing to do
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (UPstream::debug && errNo == 0)
        {
            Perr<< "UPstream::shutdown : was already finalized\n";
        }
        ourMpi = false;
        return;
    }

    if (!ourMpi)
    {
        WarningInFunction
            << "Finalizing MPI, but was initialized elsewhere\n";
    }
    ourMpi = false;


    // Abort - stop now, without any final synchonization steps!
    // -----

    if (errNo != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, errNo);
        return;
    }


    // Regular cleanup
    // ---------------

    if (UPstream::debug)
    {
        Perr<< "UPstream::shutdown\n";
    }

    // Check for any outstanding requests
    {
        label nOutstanding = 0;

        for (MPI_Request request : PstreamGlobals::outstandingRequests_)
        {
            if (MPI_REQUEST_NULL != request)
            {
                // TBD: MPI_Cancel(&request); MPI_Request_free(&request);
                ++nOutstanding;
            }
        }

        if (nOutstanding)
        {
            WarningInFunction
                << "Still have " << nOutstanding
                << " outstanding MPI requests."
                << " Should not happen for a normal code exit."
                << endl;
        }

        PstreamGlobals::outstandingRequests_.clear();
    }


    {
        detachOurBuffers();

        forAllReverse(myProcNo_, communicator)
        {
            freeCommunicatorComponents(communicator);
        }
    }


    MPI_Finalize();
}


void Foam::UPstream::exit(int errNo)
{
    UPstream::shutdown(errNo);
    std::exit(errNo);
}


void Foam::UPstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::UPstream::allocateCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{
    PstreamGlobals::initCommunicator(index);

    int returnCode = MPI_SUCCESS;

    if (parentIndex == -1)
    {
        // Global communicator. Same as world communicator for single-world

        if (index != UPstream::commGlobal())
        {
            FatalErrorInFunction
                << "base world communicator should always be index "
                << UPstream::commGlobal()
                << Foam::exit(FatalError);
        }
        auto& mpiNewComm = PstreamGlobals::MPICommunicators_[index];

        // PstreamGlobals::pendingMPIFree_[index] = false;
        // PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;

        PstreamGlobals::pendingMPIFree_[index] = true;
        MPI_Comm_dup(MPI_COMM_WORLD, &mpiNewComm);

        MPI_Comm_rank(mpiNewComm, &myProcNo_[index]);

        // Set the number of ranks to the actual number
        int numProcs;
        MPI_Comm_size(mpiNewComm, &numProcs);

        // identity [0-numProcs], as 'int'
        procIDs_[index].resize_nocopy(numProcs);
        std::iota(procIDs_[index].begin(), procIDs_[index].end(), 0);
    }
    else if (parentIndex == -2)
    {
        // MPI_COMM_SELF

        PstreamGlobals::pendingMPIFree_[index] = false;
        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_SELF;

        MPI_Comm_rank(MPI_COMM_SELF, &myProcNo_[index]);

        // For MPI_COMM_SELF : the process IDs within the world communicator.
        // Uses MPI_COMM_WORLD in case called before UPstream::commGlobal()
        // was initialized

        procIDs_[index].resize_nocopy(1);
        MPI_Comm_rank(MPI_COMM_WORLD, &procIDs_[index].front());
    }
    else
    {
        // General sub-communicator.
        // Create based on the groupings predefined by procIDs_

        const auto mpiParentComm =
            PstreamGlobals::MPICommunicators_[parentIndex];

        auto& mpiNewComm =
            PstreamGlobals::MPICommunicators_[index];

        PstreamGlobals::pendingMPIFree_[index] = true;

        // Starting from parent
        MPI_Group parent_group;
        MPI_Comm_group(mpiParentComm, &parent_group);

        MPI_Group active_group;
        MPI_Group_incl
        (
            parent_group,
            procIDs_[index].size(),
            procIDs_[index].cdata(),
           &active_group
        );

        #if defined(MSMPI_VER)
        // ms-mpi (10.0 and others?) does not have MPI_Comm_create_group
        MPI_Comm_create
        (
            mpiParentComm,
            active_group,
           &mpiNewComm
        );
        #else
        // Create new communicator for this group
        MPI_Comm_create_group
        (
            mpiParentComm,
            active_group,
            UPstream::msgType(),
           &mpiNewComm
        );
        #endif

        // Groups not needed after this...
        MPI_Group_free(&parent_group);
        MPI_Group_free(&active_group);

        if (MPI_COMM_NULL == mpiNewComm)
        {
            // This process is not involved in the new communication pattern
            myProcNo_[index] = -1;
            procIDs_[index].clear();
            PstreamGlobals::pendingMPIFree_[index] = false;
        }
        else
        {
            returnCode = MPI_Comm_rank(mpiNewComm, &myProcNo_[index]);

            if (FOAM_UNLIKELY(MPI_SUCCESS != returnCode))
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << flatOutput(procIDs_[index])
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::dupCommunicatorComponents
(
    const label parentIndex,
    const label index
)
{
    PstreamGlobals::initCommunicator(index);

    PstreamGlobals::pendingMPIFree_[index] = true;
    MPI_Comm_dup
    (
        PstreamGlobals::MPICommunicators_[parentIndex],
       &PstreamGlobals::MPICommunicators_[index]
    );

    myProcNo_[index] = myProcNo_[parentIndex];
    procIDs_[index] = procIDs_[parentIndex];
}


void Foam::UPstream::splitCommunicatorComponents
(
    const label parentIndex,
    const label index,
    int colour
)
{
    PstreamGlobals::initCommunicator(index);

    // ------------------------------------------------------------------------
    // Create sub-communicator according to its colouring
    //     => MPI_Comm_split().
    // Since other parts of OpenFOAM may still need a view of the siblings:
    //     => MPI_Group_translate_ranks().
    //
    // The MPI_Group_translate_ranks() step can be replaced with an
    // MPI_Allgather() of the involved parent ranks (since we alway maintain
    // the relative rank order when splitting).
    //
    // Since MPI_Comm_split() already does an MPI_Allgather() internally
    // to pick out the colours (and do any sorting), we can simply to
    // do the same thing:
    //
    // Do the Allgather first and pickout identical colours to define the
    // group and create a communicator based on that.
    //
    // This is no worse than the Allgather communication overhead of using
    // MPI_Comm_split() and saves the extra translate_ranks step.
    // ------------------------------------------------------------------------

    const auto mpiParentComm = PstreamGlobals::MPICommunicators_[parentIndex];

    int parentRank = 0;
    int parentSize = 0;
    MPI_Comm_rank(mpiParentComm, &parentRank);
    MPI_Comm_size(mpiParentComm, &parentSize);

    // Initialize, first marking the 'procIDs_' with the colours
    auto& procIds = procIDs_[index];

    myProcNo_[index] = -1;
    procIds.resize_nocopy(parentSize);
    procIds[parentRank] = colour;

    MPI_Allgather
    (
        MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
        procIds.data(), 1, MPI_INT,
        mpiParentComm
    );


    if (colour < 0)
    {
        procIds.clear();
    }
    else
    {
        auto last =
            std::copy_if
            (
                procIds.cbegin(),
                procIds.cend(),
                procIds.begin(),
                [=](int c){ return (c == colour); }
            );

        procIds.resize(std::distance(procIds.begin(), last));
    }

    allocateCommunicatorComponents(parentIndex, index);
}


void Foam::UPstream::freeCommunicatorComponents(const label index)
{
    if (UPstream::debug)
    {
        Perr<< "freeCommunicatorComponents: " << index
            << " from " << PstreamGlobals::MPICommunicators_.size() << endl;
    }

    // Only free communicators that we have specifically allocated ourselves
    //
    // Bounds checking needed since there are no UPstream communicator indices
    // when MPI is initialized outside of OpenFOAM

    if
    (
        (index >= 0 && index < PstreamGlobals::MPICommunicators_.size())
     && PstreamGlobals::pendingMPIFree_[index]
    )
    {
        PstreamGlobals::pendingMPIFree_[index] = false;

        // Free communicator. Sets communicator to MPI_COMM_NULL
        if (MPI_COMM_NULL != PstreamGlobals::MPICommunicators_[index])
        {
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[index]);
        }
    }
}


bool Foam::UPstream::setSharedMemoryCommunicators()
{
    // Uses the world communicator (not global communicator)

    // Skip if non-parallel
    if (!UPstream::parRun())
    {
        numNodes_ = 1;
        return false;
    }

    if (FOAM_UNLIKELY(interHostComm_ >= 0 || intraHostComm_ >= 0))
    {
        // Failed sanity check
        FatalErrorInFunction
            << "Host communicator(s) already created!" << endl
            << Foam::abort(FatalError);
        return false;
    }

    interHostComm_ = getAvailableCommIndex(UPstream::constWorldComm_);
    intraHostComm_ = getAvailableCommIndex(UPstream::constWorldComm_);

    PstreamGlobals::initCommunicator(interHostComm_);
    PstreamGlobals::initCommunicator(intraHostComm_);

    // Overwritten later
    myProcNo_[interHostComm_] = UPstream::masterNo();
    myProcNo_[intraHostComm_] = UPstream::masterNo();

    // Sorted order, purely cosmetic
    if (intraHostComm_ < interHostComm_)
    {
        std::swap(intraHostComm_, interHostComm_);
    }

    if (debug)
    {
        Perr<< "Allocating host communicators "
            << interHostComm_ << ", " << intraHostComm_ << nl
            << "    parent : " << UPstream::constWorldComm_ << nl
            << endl;
    }


    const auto mpiParentComm =
        PstreamGlobals::MPICommunicators_[UPstream::constWorldComm_];

    auto& mpiIntraNode =
        PstreamGlobals::MPICommunicators_[intraHostComm_];

    int parentRank = 0;
    int parentSize = 0;
    MPI_Comm_rank(mpiParentComm, &parentRank);
    MPI_Comm_size(mpiParentComm, &parentSize);

    List<int> hostLeaders(parentSize);
    hostLeaders = -1;

    MPI_Comm_split_type
    (
        mpiParentComm,
        MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
       &mpiIntraNode
    );

    if (FOAM_UNLIKELY(MPI_COMM_NULL == mpiIntraNode))
    {
        // This process is not involved in an intra-host communication?
        // - should never happen!

        const label index = intraHostComm_;
        PstreamGlobals::pendingMPIFree_[index] = false;

        myProcNo_[index] = -1;
        procIDs_[index].clear();

        FatalErrorInFunction
            << "Comm_split_type(shared) failed\n"
            << Foam::abort(FatalError);
    }
    else
    {
        // This process is involved in intra-host communication
        const label index = intraHostComm_;
        auto& procIds = procIDs_[index];

        PstreamGlobals::pendingMPIFree_[index] = true;

        int localRank = 0;
        int localSize = 0;
        MPI_Comm_rank(mpiIntraNode, &localRank);
        MPI_Comm_size(mpiIntraNode, &localSize);

        if (localRank == 0)
        {
            // This process is a host leader - mark its position
            hostLeaders[parentRank] = parentRank;
        }

        procIds.resize_nocopy(localSize);
        procIds[localRank] = UPstream::myProcNo(UPstream::constWorldComm_);
        // OR: procIds[localRank] = parentRank;

        // Get all of the siblings (within the node)
        MPI_Allgather
        (
            MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            procIds.data(), 1, MPI_INT,
            mpiIntraNode
        );
    }


    // Get all of the host-leader information and find who they are.
    {
        auto& procIds = procIDs_[interHostComm_];

        MPI_Allgather
        (
            MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            hostLeaders.data(), 1, MPI_INT,
            mpiParentComm
        );

        // Capture the size (number of hosts) before doing anything further
        numNodes_ = std::count_if
        (
            hostLeaders.cbegin(),
            hostLeaders.cend(),
            [](int rank){ return (rank >= 0); }
        );

        if (hostLeaders[parentRank] < 0)
        {
            // This process is not a host leader
            procIds.clear();
        }
        else
        {
            // This process is a host leader
            procIds.resize_nocopy(numNodes_);

            std::copy_if
            (
                hostLeaders.cbegin(),
                hostLeaders.cend(),
                procIds.begin(),
                [](int rank){ return (rank >= 0); }
            );
        }
    }

    // From master to host-leader. Ranks between hosts.
    allocateCommunicatorComponents(UPstream::worldComm, interHostComm_);

    return true;
}


void Foam::UPstream::barrier(const label communicator, UPstream::Request* req)
{
    // No-op for non-parallel or not on communicator
    if (!UPstream::parRun() || !UPstream::is_rank(communicator))
    {
        PstreamGlobals::reset_request(req);
        return;
    }

    if (req)
    {
        MPI_Request request;

        // Non-blocking
        if
        (
            MPI_Ibarrier
            (
                PstreamGlobals::MPICommunicators_[communicator],
               &request
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Ibarrier returned with error"
                << Foam::abort(FatalError);
        }

        *req = UPstream::Request(request);
    }
    else
    {
        // Blocking
        if
        (
            MPI_Barrier
            (
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Barrier returned with error"
                << Foam::abort(FatalError);
        }
    }
}


std::pair<int,int64_t>
Foam::UPstream::probeMessage
(
    const UPstream::commsTypes commsType,
    const int fromProcNo,
    const int tag,
    const label communicator
)
{
    std::pair<int,int64_t> result(-1, 0);

    // No-op for non-parallel or not on communicator
    if (!UPstream::parRun() || !UPstream::is_rank(communicator))
    {
        return result;
    }

    const int source = (fromProcNo < 0) ? MPI_ANY_SOURCE : fromProcNo;
    // Supporting MPI_ANY_TAG is not particularly useful...

    int flag = 0;
    MPI_Status status;

    if (UPstream::commsTypes::buffered == commsType)
    {
        // Blocking
        profilingPstream::beginTiming();

        if
        (
            MPI_Probe
            (
                source,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &status
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Probe returned with error"
                << Foam::abort(FatalError);
        }

        profilingPstream::addProbeTime();
        flag = 1;
    }
    else
    {
        // Non-blocking
        profilingPstream::beginTiming();

        if
        (
            MPI_Iprobe
            (
                source,
                tag,
                PstreamGlobals::MPICommunicators_[communicator],
                &flag,
                &status
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Iprobe returned with error"
                << Foam::abort(FatalError);
        }

        profilingPstream::addRequestTime();
    }

    if (flag)
    {
        // Unlikely to be used with large amounts of data,
        // but use MPI_Get_elements_x() instead of MPI_Count() anyhow

        #ifdef Pstream_use_MPI_Get_count
        int count(0);
        MPI_Get_count(&status, MPI_BYTE, &count);
        #else
        MPI_Count count(0);
        MPI_Get_elements_x(&status, MPI_BYTE, &count);
        #endif

        // Errors
        if (count == MPI_UNDEFINED || int64_t(count) < 0)
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "returned undefined or negative value"
                << Foam::abort(FatalError);
        }
        else if (int64_t(count) > int64_t(INT_MAX))
        {
            FatalErrorInFunction
                << "MPI_Get_count() or MPI_Get_elements_x() : "
                   "count is larger than INI_MAX bytes"
                << Foam::abort(FatalError);
        }


        result.first = status.MPI_SOURCE;
        result.second = int64_t(count);
    }

    return result;
}


// ************************************************************************* //
