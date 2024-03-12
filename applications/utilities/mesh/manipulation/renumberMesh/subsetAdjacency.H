/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
   This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Subsetting of an adjacency matrix (as CompactListList).
    Can be relocated elsewhere.

\*---------------------------------------------------------------------------*/

#include "CompactListList.H"
#include "bitSet.H"

namespace Foam
{

// Perform a subset of the adjacency matrix
CompactListList<label> subsetAdjacency
(
    const bitSet& select,
    const CompactListList<label>& input,
    labelList& subMap
)
{
    // Corresponds to cellMap etc (the original selection)
    subMap = select.toc();

    // Ensure that the subMap corresponds to a valid subset
    {
        label validSize = 0;

        const label nTotal = input.size();

        forAllReverse(subMap, i)
        {
            if (subMap[i] < nTotal)
            {
                validSize = i + 1;
                break;
            }
        }

        subMap.resize(validSize);
    }


    // Pass 1: determine the selected sub-sizes
    labelList sizes(subMap.size(), Foam::zero{});

    forAll(subMap, idx)
    {
        const auto& neighbours = input[subMap[idx]];

        for (const label nbr : neighbours)
        {
            if (select.test(nbr))
            {
                ++sizes[idx];
            }
        }
    }


    CompactListList<label> output(sizes);

    // Reuse sizes as output offset into output.values()
    sizes = labelList::subList(output.offsets(), output.size());
    labelList& values = output.values();


    // Pass 2: extract sub-adjacent matrix
    forAll(subMap, idx)
    {
        const auto& neighbours = input[subMap[idx]];

        for (const label nbr : neighbours)
        {
            if (select.test(nbr))
            {
                values[sizes[idx]++] = nbr;
            }
        }
    }

    return output;
}

} // End namespace Foam


// ************************************************************************* //
