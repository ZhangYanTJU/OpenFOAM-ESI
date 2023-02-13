/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "PstreamBuffers.H"
#ifdef Foam_PstreamBuffers_dense
#include "bitSet.H"
#endif

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef Foam_PstreamBuffers_map_storage
//- Retrieve size of specified buffer, first checking for existence
static inline label getBufferSize
(
    const Map<DynamicList<char>>& buffers,
    const label proci
)
{
    const auto iter = buffers.cfind(proci);

    return (iter.good() ? iter.val().size() : 0);
}
#else
//- Retrieve size of specified buffer, no access checking
static inline label getBufferSize
(
    const UList<DynamicList<char>>& buffers,
    const label proci
)
{
    return buffers[proci].size();
}
#endif

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PstreamBuffers::finalExchange
(
    const bool wait,
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes
    #else
    labelList& recvSizes
    #endif
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        // Use PEX algorithm

        #ifdef Foam_PstreamBuffers_map_storage
        // PEX stage 1: exchange sizes (non-blocking consensus)
        Pstream::exchangeSizes
        (
            sendBuffers_,
            recvSizes,
            tag_,
            comm_
        );
        #else
        // Like Pstream::exchangeSizes
        labelList sendSizes(nProcs_);
        forAll(sendBuffers_, proci)
        {
            sendSizes[proci] = sendBuffers_[proci].size();
        }
        recvSizes.resize_nocopy(nProcs_);

        // PEX stage 1: exchange sizes (non-blocking consensus)
        UPstream::allToAllConsensus
        (
            sendSizes,
            recvSizes,
            (tag_ + 314159),  // some unique tag?
            comm_
        );
        #endif

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


#ifdef Foam_PstreamBuffers_dense
void Foam::PstreamBuffers::finalExchange
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const bool wait,
    labelList& recvSizes
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        Pstream::exchangeSizes
        (
            sendProcs,
            recvProcs,
            sendBuffers_,
            recvSizes,
            tag_,
            comm_
        );

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}
#endif


void Foam::PstreamBuffers::finalGatherScatter
(
    const bool isGather,
    const bool wait,
    const bool needSizes,
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes
    #else
    labelList& recvSizes
    #endif
)
{
    // Could also check that it is not called twice
    // but that is used for overlapping send/recv (eg, overset)
    finishedSendsCalled_ = true;

    if (isGather)
    {
        // gather mode (all-to-one)

        // Only send to master [0]. Master is also allowed to 'send' to itself

        #ifdef Foam_PstreamBuffers_map_storage
        forAllIters(sendBuffers_, iter)
        {
            if (iter.key() != 0)
            {
                iter.val().clear();
            }
        }
        #else
        for (label proci=1; proci < sendBuffers_.size(); ++proci)
        {
            sendBuffers_[proci].clear();
        }
        #endif
    }
    else
    {
        // scatter mode (one-to-all)

        if (!UPstream::master(comm_))
        {
            // Non-master: has no sends
            clearSends();
        }
    }


    if (commsType_ == UPstream::commsTypes::nonBlocking)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        labelList recvCount;
        #else
        labelList& recvCount = recvSizes;
        #endif

        if (isGather)
        {
            // gather mode (all-to-one): master [0] <- everyone

            const label nSend = getBufferSize(sendBuffers_, 0);

            recvCount = UPstream::listGatherValues(nSend, comm_);

            #ifdef Foam_PstreamBuffers_map_storage
            // Transcribe recv count from list to map
            recvSizes.clear();
            if (UPstream::master(comm_))
            {
                for (label proci=1; proci < recvCount.size(); ++proci)
                {
                    if (recvCount[proci] > 0)
                    {
                        recvSizes.insert(proci, recvCount[proci]);
                    }
                }
            }
            #else
            if (!UPstream::master(comm_))
            {
                recvSizes.resize_nocopy(nProcs_);
                recvSizes = Zero;
            }
            #endif
        }
        else
        {
            // scatter mode (one-to-all): master [0] -> everyone

            if (UPstream::master(comm_))
            {
                recvCount.resize(nProcs_, Zero);

                #ifdef Foam_PstreamBuffers_map_storage
                forAllConstIters(sendBuffers_, iter)
                {
                    recvCount[iter.key()] = iter.val().size();
                }
                #else
                forAll(sendBuffers_, proci)
                {
                    recvCount[proci] = sendBuffers_[proci].size();
                }
                #endif
            }
            else
            {
                // Scattering, so non-master sends nothing
                recvCount = Zero;

                #ifdef Foam_PstreamBuffers_map_storage
                recvSizes.clear();
                recvSizes.resize_nocopy(nProcs_);
                #else
                recvSizes = Zero;
                #endif
            }

            const label nRecv(UPstream::listScatterValues(recvCount, comm_));

            if (UPstream::master(comm_))
            {
                #ifdef Foam_PstreamBuffers_map_storage
                recvSizes.clear();
                #else
                recvSizes = Zero;
                #endif
            }
            else
            {
                #ifdef Foam_PstreamBuffers_map_storage
                recvSizes.clear();

                if (nRecv)
                {
                    recvSizes.insert(0, nRecv);
                }
                #else
                recvSizes = Zero;
                recvSizes[0] = nRecv;
                #endif
            }
        }

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuffers_,
            recvSizes,
            recvBuffers_,
            tag_,
            comm_,
            wait
        );
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    UPstream::commsTypes commsType,
    int tag,
    label communicator,
    IOstreamOption::streamFormat fmt
)
:
    finishedSendsCalled_(false),
    allowClearRecv_(true),
    format_(fmt),
    commsType_(commsType),
    tag_(tag),
    comm_(communicator),
    nProcs_(UPstream::nProcs(comm_)),

    #ifdef Foam_PstreamBuffers_map_storage
    // Default sizing (128) is probably OK.
    // Meshes often have 16-20 neighbours (avg) and 100 neighbours (max)
    sendBuffers_(),
    recvBuffers_(),
    recvPositions_()
    #else
    sendBuffers_(nProcs_),
    recvBuffers_(nProcs_),
    recvPositions_(nProcs_, Zero)
    #endif
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::~PstreamBuffers()
{
    // Check that all data has been consumed
    #ifdef Foam_PstreamBuffers_map_storage
    forAllConstIters(recvBuffers_, iter)
    {
        const label proci = iter.key();
        const label len = iter.val().size();
        const label pos = recvPositions_.lookup(proci, len);

        if (pos < len)
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " Only consumed " << pos << " of " << len << " bytes" << nl
                << Foam::abort(FatalError);
        }
    }
    #else
    forAll(recvBuffers_, proci)
    {
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos < len)
        {
            FatalErrorInFunction
                << "Message from processor " << proci
                << " Only consumed " << pos << " of " << len << " bytes" << nl
                << Foam::abort(FatalError);
        }
    }
    #endif
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::DynamicList<char>& Foam::PstreamBuffers::accessSendBuffer
(
    const label proci
)
{
    #ifdef Foam_PstreamBuffers_map_storage
    return sendBuffers_(proci);  // Created on demand if needed
    #else
    return sendBuffers_[proci];
    #endif
}


Foam::DynamicList<char>& Foam::PstreamBuffers::accessRecvBuffer
(
    const label proci
)
{
    #ifdef Foam_PstreamBuffers_map_storage
    return recvBuffers_(proci);  // Created on demand if needed
    #else
    return recvBuffers_[proci];
    #endif
}


Foam::label& Foam::PstreamBuffers::accessRecvPosition(const label proci)
{
    #ifdef Foam_PstreamBuffers_map_storage
    return recvPositions_(proci, 0);  // Created on demand if needed
    #else
    return recvPositions_[proci];
    #endif
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::clearSends()
{
    #ifdef Foam_PstreamBuffers_map_storage
    forAllIters(sendBuffers_, iter)
    {
        iter.val().clear();
    }
    #else
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clear();
    }
    #endif
}


void Foam::PstreamBuffers::clearRecvs()
{
    #ifdef Foam_PstreamBuffers_map_storage
    forAllIters(recvBuffers_, iter)
    {
        iter.val().clear();
    }
    forAllIters(recvPositions_, iter)
    {
        iter.val() = 0;
    }
    #else
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clear();
    }
    recvPositions_ = Zero;
    #endif
}


void Foam::PstreamBuffers::clear()
{
    clearSends();
    clearRecvs();
    finishedSendsCalled_ = false;
}


void Foam::PstreamBuffers::clearSend(const label proci)
{
    #ifdef Foam_PstreamBuffers_map_storage
    {
        auto iter = sendBuffers_.find(proci);
        if (iter.good())
        {
            iter.val().clear();
        }
    }
    #else
    sendBuffers_[proci].clear();
    #endif
}


void Foam::PstreamBuffers::clearRecv(const label proci)
{
    #ifdef Foam_PstreamBuffers_map_storage
    {
        auto iter = recvBuffers_.find(proci);
        if (iter.good())
        {
            iter.val().clear();
        }
    }
    {
        auto iter = recvPositions_.find(proci);
        if (iter.good())
        {
            iter.val() = 0;
        }
    }
    #else
    recvBuffers_[proci].clear();
    recvPositions_[proci] = 0;
    #endif
}


void Foam::PstreamBuffers::clearStorage()
{
    // Could also clear out entire sendBuffers_, recvBuffers_ and reallocate.
    // Not sure if it makes much difference

    #ifdef Foam_PstreamBuffers_map_storage
    forAllIters(sendBuffers_, iter)
    {
        iter.val().clearStorage();
    }
    forAllIters(recvBuffers_, iter)
    {
        iter.val().clearStorage();
    }
    forAllIters(recvPositions_, iter)
    {
        iter.val() = 0;
    }
    #else
    for (DynamicList<char>& buf : sendBuffers_)
    {
        buf.clearStorage();
    }
    for (DynamicList<char>& buf : recvBuffers_)
    {
        buf.clearStorage();
    }
    recvPositions_ = Zero;
    #endif

    finishedSendsCalled_ = false;
}


bool Foam::PstreamBuffers::hasSendData() const
{
    #ifdef Foam_PstreamBuffers_map_storage
    forAllConstIters(sendBuffers_, iter)
    {
        if (!iter.val().empty())
        {
            return true;
        }
    }
    #else
    for (const DynamicList<char>& buf : sendBuffers_)
    {
        if (!buf.empty())
        {
            return true;
        }
    }
    #endif
    return false;
}


bool Foam::PstreamBuffers::hasRecvData() const
{
    if (finishedSendsCalled_)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        forAllConstIters(recvBuffers_, iter)
        {
            const label proci = iter.key();
            const label len = iter.val().size();

            if (recvPositions_.lookup(proci, 0) < len)
            {
                return true;
            }
        }
        #else
        forAll(recvBuffers_, proci)
        {
            if (recvPositions_[proci] < recvBuffers_[proci].size())
            {
                return true;
            }
        }
        #endif
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return false;
}


Foam::label Foam::PstreamBuffers::sendDataCount(const label proci) const
{
    return getBufferSize(sendBuffers_, proci);
}


Foam::label Foam::PstreamBuffers::recvDataCount(const label proci) const
{
    if (finishedSendsCalled_)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        #else
        const label len
        (
            getBufferSize(recvBuffers_, proci)
          - recvPositions_.lookup(proci, 0)
        );
        #else
        const label len(recvBuffers_[proci].size() - recvPositions_[proci]);
        #endif
        if (len > 0)
        {
            return len;
        }
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return 0;
}


Foam::labelList Foam::PstreamBuffers::recvDataCounts() const
{
    labelList counts(nProcs_, Zero);

    if (finishedSendsCalled_)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        forAllConstIters(recvBuffers_, iter)
        {
            const label proci = iter.key();
            const label len
            (
                iter.val().size() - recvPositions_.lookup(proci, 0)
            );

            if (len > 0)
            {
                counts[proci] = len;
            }
        }
        #else
        forAll(recvBuffers_, proci)
        {
            const label len(recvBuffers_[proci].size() - recvPositions_[proci]);

            if (len > 0)
            {
                counts[proci] = len;
            }
        }
        #endif
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return counts;
}


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount
(
    const label excludeProci
) const
{
    label maxLen = 0;

    if (finishedSendsCalled_)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        forAllConstIters(recvBuffers_, iter)
        {
            const label proci = iter.key();
            if (excludeProci != proci)
            {
                label len(iter.val().size() - recvPositions_.lookup(proci, 0));
                maxLen = max(maxLen, len);
            }
        }
        #else
        forAll(recvBuffers_, proci)
        {
            if (excludeProci != proci)
            {
                label len(recvBuffers_[proci].size() - recvPositions_[proci]);
                maxLen = max(maxLen, len);
            }
        }
        #endif
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return maxLen;
}


Foam::label Foam::PstreamBuffers::maxRecvCount() const
{
    // Use out-of-range proci to avoid excluding any processor
    return maxNonLocalRecvCount(-1);
}


Foam::label Foam::PstreamBuffers::maxNonLocalRecvCount() const
{
    return maxNonLocalRecvCount(UPstream::myProcNo(comm_));
}


const Foam::UList<char>
Foam::PstreamBuffers::peekRecvData(const label proci) const
{
    if (finishedSendsCalled_)
    {
        #ifdef Foam_PstreamBuffers_map_storage
        const auto iter = recvBuffers_.cfind(proci);

        if (iter.good())
        {
            const label pos = recvPositions_.lookup(proci, 0);
            const label len = iter.val().size();

            if (pos < len)
            {
                return UList<char>
                (
                    const_cast<char*>(iter.val().cbegin(pos)),
                    (len - pos)
                );
            }
        }
        #else
        const label pos = recvPositions_[proci];
        const label len = recvBuffers_[proci].size();

        if (pos < len)
        {
            return UList<char>
            (
                const_cast<char*>(recvBuffers_[proci].cbegin(pos)),
                (len - pos)
            );
        }
        #endif
    }
    #ifdef FULLDEBUG
    else
    {
        FatalErrorInFunction
            << "Call finishedSends first" << exit(FatalError);
    }
    #endif

    return UList<char>();
}


bool Foam::PstreamBuffers::allowClearRecv(bool on) noexcept
{
    bool old(allowClearRecv_);
    allowClearRecv_ = on;
    return old;
}


void Foam::PstreamBuffers::finishedSends(const bool wait)
{
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label> recvSizes;
    #else
    labelList recvSizes;
    #endif
    finalExchange(wait, recvSizes);
}


void Foam::PstreamBuffers::finishedSends
(
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes,
    #else
    labelList& recvSizes,
    #endif
    const bool wait
)
{
    finalExchange(wait, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


#ifdef Foam_PstreamBuffers_dense
void Foam::PstreamBuffers::finishedSends
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const bool wait
)
{
    labelList recvSizes;
    finalExchange(sendProcs, recvProcs, wait, recvSizes);
}


void Foam::PstreamBuffers::finishedSends
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    labelList& recvSizes,
    const bool wait
)
{
    finalExchange(sendProcs, recvProcs, wait, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


bool Foam::PstreamBuffers::finishedSends
(
    bitSet& sendConnections,
    DynamicList<label>& sendProcs,
    DynamicList<label>& recvProcs,
    const bool wait
)
{
    bool changed = (sendConnections.size() != nProcs());

    if (changed)
    {
        sendConnections.resize(nProcs());
    }

    // Update send connections
    // - reasonable to assume there are no self-sends on UPstream::myProcNo
    forAll(sendBuffers_, proci)
    {
        // ie, sendDataCount(proci) != 0
        if (sendConnections.set(proci, !sendBuffers_[proci].empty()))
        {
            // The state changed
            changed = true;
        }
    }

    UPstream::reduceOr(changed, comm_);

    if (changed)
    {
        // Create send/recv topology

        // The send ranks
        sendProcs.clear();
        forAll(sendBuffers_, proci)
        {
            // ie, sendDataCount(proci) != 0
            if (!sendBuffers_[proci].empty())
            {
                sendProcs.push_back(proci);
            }
        }

        labelList recvSizes;
        finishedSends(recvSizes, wait);  // All-to-all

        // The recv ranks
        recvProcs.clear();
        forAll(recvSizes, proci)
        {
            if (recvSizes[proci] > 0)
            {
                recvProcs.push_back(proci);
            }
        }
    }
    else
    {
        // Use existing send/recv ranks

        finishedSends(sendProcs, recvProcs, wait);
    }

    return changed;
}
#endif


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes,
    #else
    labelList& recvSizes,
    #endif
    const bool wait
)
{
    #ifdef Foam_PstreamBuffers_map_storage
    recvSizes.clear();

    for (const label proci : neighProcs)
    {
        recvSizes.insert(proci, 0);
    }

    // Prune any send buffers that are not neighbours
    forAllIters(sendBuffers_, iter)
    {
        if (!recvSizes.contains(iter.key()))
        {
            iter.val().clear();
        }
    }

    finalExchange(wait, recvSizes);
    #else
    // Resize for copying back
    recvSizes.resize_nocopy(sendBuffers_.size());

    // Prune send buffers that are not neighbours
    {
        labelHashSet keepProcs(neighProcs);

        // Prune send buffers that are not neighbours
        forAll(sendBuffers_, proci)
        {
            if (!keepProcs.contains(proci))
            {
                sendBuffers_[proci].clear();
            }
        }
    }

    finalExchange(wait, recvSizes);
    #endif
}


void Foam::PstreamBuffers::finishedNeighbourSends
(
    const labelUList& neighProcs,
    const bool wait
)
{
    #ifdef Foam_PstreamBuffers_map_storage
    finishedSends(neighProcs, neighProcs, wait);
    #else
    labelList recvSizes;

    // Prune send buffers that are not neighbours
    {
        labelHashSet keepProcs(neighProcs);

        // Prune send buffers that are not neighbours
        forAll(sendBuffers_, proci)
        {
            if (!keepProcs.contains(proci))
            {
                sendBuffers_[proci].clear();
            }
        }
    }

    finalExchange(wait, recvSizes);
    #endif
}


void Foam::PstreamBuffers::finishedGathers(const bool wait)
{
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label> recvSizes;
    finalGatherScatter(true, wait, false, recvSizes);
    #else
    labelList recvSizes;
    finalGatherScatter(true, wait, false, recvSizes);
    #endif
}


void Foam::PstreamBuffers::finishedScatters(const bool wait)
{
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label> recvSizes;
    finalGatherScatter(false, wait, false, recvSizes);
    #else
    labelList recvSizes;
    finalGatherScatter(false, wait, false, recvSizes);
    #endif
}


void Foam::PstreamBuffers::finishedGathers
(
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes,
    #else
    labelList& recvSizes,
    #endif
    const bool wait
)
{
    finalGatherScatter(true, wait, true, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


void Foam::PstreamBuffers::finishedScatters
(
    #ifdef Foam_PstreamBuffers_map_storage
    Map<label>& recvSizes,
    #else
    labelList& recvSizes,
    #endif
    const bool wait
)
{
    finalGatherScatter(false, wait, true, recvSizes);

    if (commsType_ != UPstream::commsTypes::nonBlocking)
    {
        FatalErrorInFunction
            << "Obtaining sizes not supported in "
            << UPstream::commsTypeNames[commsType_] << endl
            << " since transfers already in progress. Use non-blocking instead."
            << exit(FatalError);

        // Note: maybe possible only if using different tag from write started
        // by ~UOPstream. Needs some work.
    }
}


// ************************************************************************* //
