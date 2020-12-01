/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "distanceSurface.H"
#include "regionSplit.H"
#include "syncTools.H"
#include "tetMatcher.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::distanceSurface::refineBlockedCells
(
    bitSet& ignoreCells,
    const isoSurfaceBase& isoCutter
) const
{
    // With the cell/point distance fields we can take a second pass at
    // pre-filtering.
    // This duplicates how cut detection is determined in the cell/topo
    // algorithms but is fairly inexpensive (creates no geometry)

    bool changed = false;

    tetMatcher tet;

    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        if (ignoreCells.test(celli))
        {
            continue;
        }

        auto cut = isoCutter.getCellCutType(celli, tet.isA(mesh_, celli));
        if (!(cut & isoSurfaceBase::ANYCUT))
        {
            ignoreCells.set(celli);
            changed = true;
        }
    }

    return changed;
}


Foam::boolList Foam::distanceSurface::filterPrepareRegionSplit
(
    const bitSet& ignoreCells
) const
{
    // Prepare for region split

    boolList blockedFaces(mesh_.nFaces(), false);

    const labelList& faceOwn = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();

    // Could be more efficient
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        // If only one cell is blocked, the face corresponds
        // to an exposed subMesh face

        if
        (
            ignoreCells.test(faceOwn[facei])
         != ignoreCells.test(faceNei[facei])
        )
        {
            blockedFaces.set(facei);
        }
    }

    for (const polyPatch& patch : mesh_.boundaryMesh())
    {
        if (!patch.coupled())
        {
            continue;
        }

        forAll(patch, patchFacei)
        {
            const label facei = patch.start() + patchFacei;
            if (ignoreCells.test(faceOwn[facei]))
            {
                blockedFaces.set(facei);
            }
        }
    }

    syncTools::syncFaceList(mesh_, blockedFaces, xorEqOp<bool>());

    return blockedFaces;
}


void Foam::distanceSurface::filterKeepLargestRegion
(
    bitSet& ignoreCells
) const
{
    // For region split
    boolList blockedFaces
    (
        filterPrepareRegionSplit(ignoreCells)
    );


    // Split region
    regionSplit rs(mesh_, blockedFaces);
    blockedFaces.clear();

    const labelList& regionColour = rs;

    // Identical number of regions on all processors
    labelList nCutsPerRegion(rs.nRegions(), Zero);

    // Count cut cells (ie, unblocked)
    forAll(regionColour, celli)
    {
        if (!ignoreCells.test(celli))
        {
            ++nCutsPerRegion[regionColour[celli]];
        }
    }

    // Sum totals from all processors
    Pstream::listCombineGather(nCutsPerRegion, plusEqOp<label>());


    // Define which regions to keep
    boolList keepRegion(rs.nRegions(), false);

    if (Pstream::master())
    {
        const label largest = findMax(nCutsPerRegion);

        if (largest == -1)
        {
            // Should not happen
            keepRegion = true;
        }
        else
        {
            keepRegion[largest] = true;
        }

        if (debug)
        {
            Info<< "Had " << sum(nCutsPerRegion) << " cuts, in "
                << nCutsPerRegion.size() << " regions, largest is "
                << largest <<  ": " << flatOutput(nCutsPerRegion) << nl;
        }
    }

    // Send to everyone
    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            for (const int proci : Pstream::subProcs())
            {
                OPstream send(Pstream::commsTypes::scheduled, proci);
                send << keepRegion;
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            recv >> keepRegion;
        }
    }


    // Apply
    forAll(regionColour, celli)
    {
        if (!keepRegion.test(regionColour[celli]))
        {
            ignoreCells.set(celli);
        }
    }
}


void Foam::distanceSurface::filterKeepNearestRegions
(
    bitSet& ignoreCells
) const
{
    if (nearestPoints_.empty())
    {
        WarningInFunction
            << "Ignoring nearestPoints - no points provided" << nl
            << endl;
        return;
    }

    // For region split
    boolList blockedFaces
    (
        filterPrepareRegionSplit(ignoreCells)
    );

    // Split region
    regionSplit rs(mesh_, blockedFaces);
    blockedFaces.clear();

    const labelList& regionColour = rs;

    const pointField& cc = mesh_.cellCentres();
    const pointField& nearPts = nearestPoints_;

    // The magSqr distance and region
    typedef Tuple2<scalar, label> nearInfo;
    List<nearInfo> nearest(nearPts.size(), nearInfo(GREAT, -1));

    // Also collect cuts per region, may be useful for rejecting
    // small regions. Code as per filterKeepLargestRegion
    labelList nCutsPerRegion(rs.nRegions(), Zero);

    forAll(cc, celli)
    {
        if (ignoreCells.test(celli))
        {
            continue;
        }

        const point& pt = cc[celli];
        const label regioni = regionColour[celli];

        ++nCutsPerRegion[regioni];

        label pointi = 0;
        for (nearInfo& near : nearest)
        {
            const scalar distSqr = magSqr(nearPts[pointi] - pt);
            ++pointi;

            if (distSqr < near.first())
            {
                near.first() = distSqr;
                near.second() = regioni;
            }
        }
    }

    // Sum totals from all processors
    Pstream::listCombineGather(nCutsPerRegion, plusEqOp<label>());

    // Get nearest
    Pstream::listCombineGather(nearest, minFirstEqOp<scalar>());


    // Define which regions to keep

    boolList keepRegion(rs.nRegions(), false);

    if (Pstream::master())
    {
        const label largest = findMax(nCutsPerRegion);

        for (const nearInfo& near : nearest)
        {
            const scalar distSqr = near.first();
            const label regioni = near.second();

            if (regioni != -1 && distSqr < maxDistanceSqr_)
            {
                keepRegion[regioni] = true;
            }
        }

        if (debug)
        {
            Info<< "Had " << sum(nCutsPerRegion) << " cuts, in "
                << nCutsPerRegion.size() << " regions, largest is "
                << largest <<  ": " << flatOutput(nCutsPerRegion) << nl;

            Info<< "nearestPoints (max distance = "
                << sqrt(maxDistanceSqr_) << ")" << nl;

            forAll(nearPts, pointi)
            {
                const scalar dist = sqrt(nearest[pointi].first());
                const label regioni = nearest[pointi].second();

                Info<< "    " << nearPts[pointi] << " region "
                    << regioni << " distance "
                    << dist;

                if (!keepRegion.test(regioni))
                {
                    Info<< " too far";
                }
                Info<< nl;
            }
        }
    }

    // Send to everyone
    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            for (const int proci : Pstream::subProcs())
            {
                OPstream send(Pstream::commsTypes::scheduled, proci);
                send << keepRegion;
            }
        }
        else
        {
            IPstream recv
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            recv >> keepRegion;
        }
    }

    // Apply
    forAll(regionColour, celli)
    {
        if (!keepRegion.test(regionColour[celli]))
        {
            ignoreCells.set(celli);
        }
    }
}


// ************************************************************************* //
