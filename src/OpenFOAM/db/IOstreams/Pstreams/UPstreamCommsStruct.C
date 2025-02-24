/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "UPstream.H"

#include <algorithm>
#include <numeric>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// This outputs as depth-first, but graphviz sorts that for us
static void printGraph_impl
(
    Ostream& os,
    const UPstream::commsStructList& comms,
    const int proci,
    int depth,
    const int maxDepth = 1024
)
{
    if (proci >= comms.size())
    {
        // Corner case when only a single rank involved
        // (eg, for node-local communicator)
        return;
    }

    const auto& below = comms[proci].below();

    if (proci == 0)
    {
        os << nl << "// communication graph:" << nl;
        os.beginBlock("graph");

        // Prefer left-to-right layout for large graphs
        os << indent << "rankdir=LR" << nl;
    }


    // Output the immediate neighbours below

    if (below.empty())
    {
        if (proci == 0)
        {
            // A graph with a single-node (eg, self-comm)
            os << indent << proci << nl;
        }
    }
    else
    {
        os << indent << proci << " -- " << token::BEGIN_BLOCK;

        // Accumulate into ranges whenever possible
        IntRange<int> range;

        // Print accumulated range and reset
        auto emit_range = [&]()
        {
            if (!range.empty())
            {
                os << ' ';
                if (range.min() < range.max())
                {
                    os << '"' << range.min() << ".." << range.max() << '"';
                }
                else
                {
                    os << range.min();
                }
                range.reset();
            }
        };

        for (const auto nbrProci : below)
        {
            const bool terminal = comms[nbrProci].below().empty();

            if
            (
                terminal
             && (!range.empty() && (range.max()+1 == nbrProci))
            )
            {
                // Accumulate
                ++range;
                continue;
            }

            // Emit accumulated range
            emit_range();

            if (terminal)
            {
                range.reset(nbrProci, 1);
            }
            else
            {
                os << token::SPACE << nbrProci;
            }
        }

        // Emit accumulated range
        emit_range();

        os << token::SPACE << token::END_BLOCK << nl;
    }


    // Recurse into below neighbours, but limit the maximum depth
    ++depth;
    if (depth >= maxDepth && (proci != 0))
    {
        return;
    }

    for (const auto nbrProci : below)
    {
        // if (proci == nbrProci) continue;  // Extreme safety!
        printGraph_impl(os, comms, nbrProci, depth, maxDepth);
    }

    if (proci == 0)
    {
        os.endBlock();
        os << "// end graph" << nl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create a tree-like schedule. For 8 procs:
// (level 0)
//      0 receives from 1
//      2 receives from 3
//      4 receives from 5
//      6 receives from 7
// (level 1)
//      0 receives from 2
//      4 receives from 6
// (level 2)
//      0 receives from 4
//
// The sends/receives for all levels are collected per processor
// (one send per processor; multiple receives possible) creating
// a table:
//
// So per processor:
// proc     receives from   sends to
// ----     -------------   --------
//  0       1,2,4           -
//  1       -               0
//  2       3               0
//  3       -               2
//  4       5               0
//  5       -               4
//  6       7               4
//  7       -               6

namespace Foam
{

static int simpleTree
(
    const int myProci,
    const int numProcs,

    DynamicList<int>& below,
    DynamicList<int>& allBelow
)
{
    int above(-1);

    for (int mod = 2, step = 1; step < numProcs; step = mod)
    {
        mod = step * 2;

        if (myProci % mod)
        {
            // The rank above
            above = myProci - (myProci % mod);
            break;
        }
        else
        {
            for
            (
                int i = myProci + step;
                i < numProcs && i < myProci + mod;
                i += step
            )
            {
                below.push_back(i);
            }
            for
            (
                int i = myProci + step;
                i < numProcs && i < myProci + mod;
                ++i
            )
            {
                allBelow.push_back(i);
            }
        }
    }

    return above;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPstream::commsStruct::commsStruct
(
    const int above,
    List<int>&& below,
    List<int>&& allBelow,
    List<int>&& allNotBelow
)
:
    above_(above),
    below_(std::move(below)),
    allBelow_(std::move(allBelow)),
    allNotBelow_(std::move(allNotBelow))
{}


Foam::UPstream::commsStruct::commsStruct
(
    const int numProcs,
    const int myProcID,
    const int above,
    const UList<int>& below,
    const UList<int>& allBelow
)
:
    above_(above),
    below_(below),
    allBelow_(allBelow),
    allNotBelow_(numProcs - allBelow.size() - 1)
{
    List<bool> isNotBelow(numProcs, true);

    // Exclude self
    isNotBelow[myProcID] = false;

    // Exclude allBelow
    for (const auto proci : allBelow)
    {
        isNotBelow[proci] = false;
    }

    // Compacting to obtain allNotBelow_
    int nNotBelow = 0;
    for (int proci = 0; proci < numProcs; ++proci)
    {
        if (isNotBelow[proci])
        {
            allNotBelow_[nNotBelow++] = proci;
        }
    }

    if (nNotBelow != allNotBelow_.size())
    {
        FatalErrorInFunction
            << "Problem: " << nNotBelow << " != " << allNotBelow_.size() << nl
            << Foam::abort(FatalError);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::UPstream::commsStructList::printGraph
(
    Ostream& os,
    const int proci
) const
{
    // Print graph - starting at depth 0
    // Avoid corner case when only a single rank involved
    // (eg, for node-local communicator)

    if (proci < size())
    {
        printGraph_impl(os, *this, proci, 0);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::UPstream::commsStruct::nProcs() const noexcept
{
    return (1 + int(allBelow_.size() + allNotBelow_.size()));
}


void Foam::UPstream::commsStruct::reset()
{
    above_ = -1;
    below_.clear();
    allBelow_.clear();
    allNotBelow_.clear();
}


void Foam::UPstream::commsStruct::reset_linear
(
    const int myProci,
    const int numProcs
)
{
    reset();

    // Linear (flat) communication pattern
    int above(-1);
    List<int> below;

    if (myProci == 0)
    {
        below.resize(numProcs-1);
        std::iota(below.begin(), below.end(), 1);
    }
    else
    {
        above = 0;
    }

    *this = UPstream::commsStruct(numProcs, myProci, above, below, below);
}


void Foam::UPstream::commsStruct::reset
(
    const int myProci,
    const int numProcs,
    const int communicator
)
{
    // Trivially small domains
    if (numProcs <= 2)
    {
        reset_linear(myProci, numProcs);
        return;
    }


    reset();

    int above(-1);
    DynamicList<int> below;
    DynamicList<int> allBelow;

    if (UPstream::usingNodeComms(communicator))
    {
        // Additional treatment...
    }


    // Simple tree communication pattern
    above = simpleTree
    (
        myProci,
        numProcs,
        below,
        allBelow
    );

    *this = UPstream::commsStruct(numProcs, myProci, above, below, allBelow);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::UPstream::commsStructList&
Foam::UPstream::commsStructList::null()
{
    static std::unique_ptr<commsStructList> singleton;

    if (!singleton)
    {
        singleton = std::make_unique<commsStructList>();
    }

    return *singleton;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::UPstream::commsStructList::linear(bool on)
{
    if (flat_ != on)
    {
        // Current size
        const auto len = tree_.size();

        flat_ = on;
        tree_.clear();
        if (len > 0)
        {
            tree_.resize(len);
        }
    }
}


void Foam::UPstream::commsStructList::reset(int communicator)
{
    comm_ = communicator;
    tree_.clear();
}


void Foam::UPstream::commsStructList::reset(int communicator, bool flat)
{
    comm_ = communicator;
    tree_.clear();
    flat_ = flat;
}


void Foam::UPstream::commsStructList::init(int communicator)
{
    comm_ = communicator;
    tree_.clear();
    if (comm_ >= 0)
    {
        tree_.resize(UPstream::nProcs(comm_));
    }
}


void Foam::UPstream::commsStructList::init(int communicator, bool flat)
{
    init(communicator);
    flat_ = flat;
}


const Foam::UPstream::commsStruct&
Foam::UPstream::commsStructList::get(int proci) const
{
    const auto numProcs = UPstream::nProcs(comm_);

    // Only if reset(comm) instead of init(comm) was used
    if (tree_.size() < numProcs)
    {
        const_cast<List<commsStruct>&>(tree_).resize(numProcs);
    }

    const UPstream::commsStruct& entry = tree_[proci];

    if (entry.nProcs() != numProcs)
    {
        // Create/update

        if (flat_)
        {
            const_cast<UPstream::commsStruct&>(entry)
                .reset_linear(proci, numProcs);
        }
        else
        {
            const_cast<UPstream::commsStruct&>(entry)
                .reset(proci, numProcs, comm_);
        }
    }

    return entry;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::UPstream::commsStruct::operator==(const commsStruct& comm) const
{
    return
    (
        (above() == comm.above())
     && (below() == comm.below())
    );
}


bool Foam::UPstream::commsStruct::operator!=(const commsStruct& comm) const
{
    return !operator==(comm);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const UPstream::commsStruct& comm)
{
    os  << comm.above() << nl;
    os  << "  "; comm.below().writeList(os) << nl;
    os  << "  "; comm.allBelow().writeList(os) << nl;
    os  << "  "; comm.allNotBelow().writeList(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
