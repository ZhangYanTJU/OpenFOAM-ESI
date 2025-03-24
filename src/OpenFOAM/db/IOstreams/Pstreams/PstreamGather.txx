/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Gather data from all processors onto single processor according to some
    communication schedule (usually tree-to-master).
    The gathered data will be a single value constructed from the values
    on individual processors using a user-specified operator.

Note
    Normal gather uses:
    - binary operator that returns a value.
      So assignment that return value to yield the new value

    Combine gather uses:
    - binary operator modifies its first parameter in-place

\*---------------------------------------------------------------------------*/

#include "contiguous.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Single value variants

template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::gather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    T& value,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const label myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();


        // Receive from my downstairs neighbours
        for (const auto proci : below)
        {
            T received;

            if constexpr (is_contiguous_v<T>)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    reinterpret_cast<char*>(&received),
                    sizeof(T),
                    tag,
                    communicator
                );
            }
            else
            {
                IPstream::recv(received, proci, tag, communicator);
            }

            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " received from "
                        << proci << " data:" << received << endl;
                }
            }

            if constexpr (InplaceMode)
            {
                // In-place binary operation
                bop(value, received);
            }
            else
            {
                // Assign result of binary operation
                value = bop(value, received);
            }
        }

        // Send up value
        if (myComm.above() >= 0)
        {
            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " sending to " << myComm.above()
                        << " data:" << value << endl;
                }
            }

            if constexpr (is_contiguous_v<T>)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    reinterpret_cast<const char*>(&value),
                    sizeof(T),
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream::send(value, myComm.above(), tag, communicator);
            }
        }
    }
}


template<class T, class BinaryOp, bool InplaceMode>
bool Foam::Pstream::gather_topo_algorithm
(
    T& value,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator)
     && UPstream::usingTopoControl(UPstream::topoControls::combine)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological reduce
        // - linear for local-node (assume communication is fast)
        // - tree for inter-node (no assumption about speed)

        using control = std::pair<int, bool>;

        for
        (
            auto [subComm, linear] :
            {
                // 1: within node
                control{ UPstream::commLocalNode(), true },
                // 2: between nodes
                control{ UPstream::commInterNode(), false }
            }
        )
        {
            if (UPstream::is_parallel(subComm))
            {
                Pstream::gather_algorithm<T, BinaryOp, InplaceMode>
                (
                    UPstream::whichCommunication(subComm, linear),
                    value,
                    bop,
                    tag,
                    subComm
                );
            }
        }
    }

    return withTopo;
}


template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::gather
(
    T& value,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiReduce
        (
            &value,
            1,
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else if
    (
        !Pstream::gather_topo_algorithm<T, BinaryOp, InplaceMode>
        (
            value,
            bop,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        Pstream::gather_algorithm<T, BinaryOp, InplaceMode>
        (
            commOrder,
            value,
            bop,
            tag,
            communicator
        );
    }
}


template<class T, class CombineOp>
void Foam::Pstream::combineGather
(
    T& value,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::gather<T, CombineOp, true>(value, cop, tag, comm);
}


template<class T, class CombineOp>
void Foam::Pstream::combineReduce
(
    T& value,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    if (UPstream::is_parallel(comm))
    {
        // In-place binary operation
        Pstream::gather<T, CombineOp, true>(value, cop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// List variants

template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::listGather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    UList<T>& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator) || values.empty())
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const label myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();

        // Same length on all ranks
        const label listLen = values.size();

        List<T> received;

        if (!below.empty())
        {
            // Pre-size for contiguous reading
            if constexpr (is_contiguous_v<T>)
            {
                received.resize_nocopy(listLen);
            }
        }

        // Receive from my downstairs neighbours
        for (const auto proci : below)
        {
            if constexpr (is_contiguous_v<T>)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    received,
                    tag,
                    communicator
                );
            }
            else
            {
                received.clear();  // extra safety?
                IPstream::recv(received, proci, tag, communicator);
            }

            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " received from "
                        << proci << " data:" << received << endl;
                }
            }

            for (label i = 0; i < listLen; ++i)
            {
                if constexpr (InplaceMode)
                {
                    // In-place binary operation
                    bop(values[i], received[i]);
                }
                else
                {
                    // Assign result of binary operation
                    values[i] = bop(values[i], received[i]);
                }
            }
        }

        // Send up values
        if (myComm.above() >= 0)
        {
            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " sending to " << myComm.above()
                        << " data:" << values << endl;
                }
            }

            if constexpr (is_contiguous_v<T>)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::scheduled,
                    myComm.above(),
                    values,
                    tag,
                    communicator
                );
            }
            else
            {
                OPstream::send(values, myComm.above(), tag, communicator);
            }
        }
    }
}


template<class T, class BinaryOp, bool InplaceMode>
bool Foam::Pstream::listGather_topo_algorithm
(
    UList<T>& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator) && !values.empty()
     && UPstream::usingTopoControl(UPstream::topoControls::combine)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological reduce
        // - linear for local-node (assume communication is fast)
        // - tree for inter-node (no assumption about speed)

        using control = std::pair<int, bool>;

        for
        (
            auto [subComm, linear] :
            {
                // 1: within node
                control{ UPstream::commLocalNode(), true },
                // 2: between nodes
                control{ UPstream::commInterNode(), false }
            }
        )
        {
            if (UPstream::is_parallel(subComm))
            {
                Pstream::listGather_algorithm<T, BinaryOp, InplaceMode>
                (
                    UPstream::whichCommunication(subComm, linear),
                    values,
                    bop,
                    tag,
                    subComm
                );
            }
        }
    }

    return withTopo;
}


template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::listGather
(
    UList<T>& values,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator) || values.empty())
    {
        // Nothing to do
        return;
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiReduce
        (
            values.data(),
            values.size(),  // Same length on all ranks
            UPstream_opType<BinaryOp>::opcode_id,
            communicator
        );
    }
    else if (values.size() == 1)
    {
        // Single value - optimized version
        Pstream::gather<T, BinaryOp, InplaceMode>
        (
            values[0],
            bop,
            tag,
            communicator
        );
    }
    else if
    (
        !Pstream::listGather_topo_algorithm<T, BinaryOp, InplaceMode>
        (
            values,
            bop,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        Pstream::listGather_algorithm<T, BinaryOp, InplaceMode>
        (
            commOrder,
            values,
            bop,
            tag,
            communicator
        );
    }
}


template<class T, class BinaryOp, bool InplaceMode>
void Foam::Pstream::listReduce
(
    UList<T>& values,
    [[maybe_unused]] BinaryOp bop,
    [[maybe_unused]] const int tag,
    const int comm
)
{
    if (!UPstream::is_parallel(comm) || values.empty())
    {
        // Nothing to do
    }
    else if constexpr (!InplaceMode && UPstream_data_opType<BinaryOp, T>::value)
    {
        // Valid opcode and (directly/indirectly) uses basic dataType
        UPstream::mpiAllReduce
        (
            values.data(),
            values.size(),  // Same length on all ranks
            UPstream_opType<BinaryOp>::opcode_id,
            comm
        );
    }
    else if (values.size() == 1)
    {
        // Single value - optimized version
        Pstream::gather<T, BinaryOp, InplaceMode>(values[0], bop, tag, comm);
        Pstream::broadcast(values[0], comm);
    }
    else
    {
        // Multiple values
        Pstream::listGather<T, BinaryOp, InplaceMode>(values, bop, tag, comm);
        Pstream::broadcast(values, comm);
    }
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineGather
(
    UList<T>& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::listGather<T, CombineOp, true>(values, cop, tag, comm);
}


template<class T, class CombineOp>
void Foam::Pstream::listCombineReduce
(
    UList<T>& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::listReduce<T, CombineOp, true>(values, cop, tag, comm);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Map variants

template<class Container, class BinaryOp, bool InplaceMode>
void Foam::Pstream::mapGather_algorithm
(
    const UPstream::commsStructList& comms,  // Communication order
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else
    {
        // if (comms.empty()) return;  // extra safety?
        const label myProci = UPstream::myProcNo(communicator);
        const auto& myComm = comms[myProci];
        const auto& below = myComm.below();


        // Receive from my downstairs neighbours
        for (const auto proci : below)
        {
            // Map/HashTable: non-contiguous
            Container received;
            IPstream::recv(received, proci, tag, communicator);

            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " received from "
                        << proci << " data:" << received << endl;
                }
            }

            const auto last = received.end();

            for (auto iter = received.begin(); iter != last; ++iter)
            {
                auto slot = values.find(iter.key());

                if (slot.good())
                {
                    // Combine with existing entry

                    if constexpr (InplaceMode)
                    {
                        // In-place binary operation
                        bop(slot.val(), iter.val());
                    }
                    else
                    {
                        // Assign result of binary operation
                        slot.val() = bop(slot.val(), iter.val());
                    }
                }
                else
                {
                    // Create a new entry
                    values.emplace(iter.key(), std::move(iter.val()));
                }

            }
        }

        // Send up values
        if (myComm.above() >= 0)
        {
            if constexpr (InplaceMode)
            {
                if (debug & 2)
                {
                    Perr<< " sending to " << myComm.above()
                        << " data:" << values << endl;
                }
            }

            OPstream::send(values, myComm.above(), tag, communicator);
        }
    }
}


template<class Container, class BinaryOp, bool InplaceMode>
bool Foam::Pstream::mapGather_topo_algorithm
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    const bool withTopo =
    (
        UPstream::is_parallel(communicator)
     && UPstream::usingTopoControl(UPstream::topoControls::mapGather)
     && UPstream::usingNodeComms(communicator)
    );

    if (withTopo)
    {
        // Topological reduce
        // - linear for local-node (assume communication is fast)
        // - tree for inter-node (no assumption about speed)

        using control = std::pair<int, bool>;

        for
        (
            auto [subComm, linear] :
            {
                // 1: within node
                control{ UPstream::commLocalNode(), true },
                // 2: between nodes
                control{ UPstream::commInterNode(), false }
            }
        )
        {
            if (UPstream::is_parallel(subComm))
            {
                Pstream::mapGather_algorithm<Container, BinaryOp, InplaceMode>
                (
                    UPstream::whichCommunication(subComm, linear),
                    values,
                    bop,
                    tag,
                    subComm
                );
            }
        }
    }

    return withTopo;
}


template<class Container, class BinaryOp, bool InplaceMode>
void Foam::Pstream::mapGather
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int communicator
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // Nothing to do
        return;
    }
    else if
    (
        !Pstream::mapGather_topo_algorithm<Container, BinaryOp, InplaceMode>
        (
            values,
            bop,
            tag,
            communicator
        )
    )
    {
        // Communication order
        const auto& commOrder = UPstream::whichCommunication(communicator);

        Pstream::mapGather_algorithm<Container, BinaryOp, InplaceMode>
        (
            commOrder,
            values,
            bop,
            tag,
            communicator
        );
    }
}


template<class Container, class BinaryOp, bool InplaceMode>
void Foam::Pstream::mapReduce
(
    Container& values,
    BinaryOp bop,
    const int tag,
    const int comm
)
{
    Pstream::mapGather<Container, BinaryOp, InplaceMode>
    (
        values, bop, tag, comm
    );
    Pstream::broadcast(values, comm);
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineGather
(
    Container& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::mapGather<Container, CombineOp, true>
    (
        values, cop, tag, comm
    );
}


template<class Container, class CombineOp>
void Foam::Pstream::mapCombineReduce
(
    Container& values,
    CombineOp cop,
    const int tag,
    const int comm
)
{
    // In-place binary operation
    Pstream::mapReduce<Container, CombineOp, true>
    (
        values, cop, tag, comm
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Single values to/from a list

template<class T>
Foam::List<T> Foam::Pstream::listGatherValues
(
    const T& localValue,
    const int communicator,
    [[maybe_unused]] const int tag
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
    else if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listGatherValues(localValue, communicator);
    }
    else
    {
        // Standard gather (all to one)

        // The data are non-contiguous!
        //
        // Non-trivial to manage non-blocking gather without a
        // PEX/NBX approach (eg, PstreamBuffers).
        // Leave with simple exchange for now

        List<T> allValues;
        if (UPstream::master(communicator))
        {
            allValues.resize(UPstream::nProcs(communicator));

            for (const int proci : UPstream::subProcs(communicator))
            {
                IPstream::recv(allValues[proci], proci, tag, communicator);
            }

            allValues[0] = localValue;
        }
        else if (UPstream::is_rank(communicator))
        {
            OPstream::send(localValue, UPstream::masterNo(), tag, communicator);
        }

        return allValues;
    }
}


template<class T>
T Foam::Pstream::listScatterValues
(
    const UList<T>& allValues,
    const int communicator,
    [[maybe_unused]] const int tag
)
{
    if (!UPstream::is_parallel(communicator))
    {
        // non-parallel: return first value
        // TBD: only when UPstream::is_rank(communicator) as well?

        if (!allValues.empty())
        {
            return allValues[0];
        }

        return T{};  // Fallback value
    }
    else if constexpr (is_contiguous_v<T>)
    {
        // UPstream version is contiguous only
        return UPstream::listScatterValues(allValues, communicator);
    }
    else
    {
        // Standard scatter (one to all)

        T localValue{};

        if (UPstream::master(communicator))
        {
            const label numProc = UPstream::nProcs(communicator);

            if (allValues.size() < numProc)
            {
                FatalErrorInFunction
                    << "Attempting to send " << allValues.size()
                    << " values to " << numProc << " processors" << endl
                    << Foam::abort(FatalError);
            }

            const label startOfRequests = UPstream::nRequests();

            List<DynamicList<char>> sendBuffers(numProc);

            for (const int proci : UPstream::subProcs(communicator))
            {
                UOPstream toProc
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    sendBuffers[proci],
                    tag,
                    communicator
                );
                toProc << allValues[proci];
            }

            // Wait for outstanding requests
            UPstream::waitRequests(startOfRequests);

            return allValues[0];
        }
        else if (UPstream::is_rank(communicator))
        {
            IPstream::recv(localValue, UPstream::masterNo(), tag, communicator);
        }

        return localValue;
    }
}


// ************************************************************************* //
