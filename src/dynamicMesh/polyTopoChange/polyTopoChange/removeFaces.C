/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "removeFaces.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyRemoveCell.H"
#include "polyRemovePoint.H"
#include "syncTools.H"
#include "OFstream.H"
#include "indirectPrimitivePatch.H"
#include "Time.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removeFaces, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Changes region of connected set of cells. Can be recursive since hopefully
// only small area of faces removed in one go.
void Foam::removeFaces::changeCellRegion
(
    const label celli,
    const label oldRegion,
    const label newRegion,
    labelList& cellRegion
) const
{
    if (cellRegion[celli] == oldRegion)
    {
        cellRegion[celli] = newRegion;

        // Step to neighbouring cells

        const labelList& cCells = mesh_.cellCells()[celli];

        forAll(cCells, i)
        {
            changeCellRegion(cCells[i], oldRegion, newRegion, cellRegion);
        }
    }
}


Foam::label Foam::removeFaces::changeFaceRegion
(
    const labelList& cellRegion,
    const boolList& removedFace,
    const labelList& nFacesPerEdge,
    const label faceI,
    const label newRegion,

    labelList& faceRegion
) const
{
    // Count number of changed faces
    label nChanged = 0;

    if (faceRegion[faceI] == -1 && !removedFace[faceI])
    {
        faceRegion[faceI] = newRegion;

        nChanged = 1;

        // Get mesh data
        const labelListList& meshFaceEdges = mesh_.faceEdges();
        const labelListList& meshEdgeFaces = mesh_.edgeFaces();

        // Step to neighbouring faces across edges that will get removed
        const labelList& fEdges = meshFaceEdges[faceI];

        forAll(fEdges, i)
        {
            const label& edgeI = fEdges[i];

            if (nFacesPerEdge[edgeI] >= 0 && nFacesPerEdge[edgeI] <= 2)
            {
                const labelList& eFaces = meshEdgeFaces[edgeI];

                forAll(eFaces, j)
                {
                    nChanged += changeFaceRegion
                    (
                        cellRegion,
                        removedFace,
                        nFacesPerEdge,
                        eFaces[j],
                        newRegion,

                        faceRegion
                    );
                }
            }
        }
    }

    return nChanged;
}


// Changes region of connected set of faces. Returns number of changed faces.
Foam::label Foam::removeFaces::changeFaceRegion
(
    const labelList& cellRegion,
    const boolList& removedFace,
    const labelList& nFacesPerEdge,
    const label facei,
    const label newRegion,
    const labelList& fEdges,
    labelList& faceRegion
) const
{
    label nChanged = 0;

    if (faceRegion[facei] == -1 && !removedFace[facei])
    {
        faceRegion[facei] = newRegion;

        nChanged = 1;

        // Storage for on-the-fly addressing
        DynamicList<label> fe;
        DynamicList<label> ef;

        // Step to neighbouring faces across edges that will get removed
        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            if (nFacesPerEdge[edgeI] >= 0 && nFacesPerEdge[edgeI] <= 2)
            {
                const labelList& eFaces = mesh_.edgeFaces(edgeI, ef);

                forAll(eFaces, j)
                {
                    label nbrFacei = eFaces[j];

                    const labelList& fEdges1 = mesh_.faceEdges(nbrFacei, fe);

                    nChanged += changeFaceRegion
                    (
                        cellRegion,
                        removedFace,
                        nFacesPerEdge,
                        nbrFacei,
                        newRegion,
                        fEdges1,
                        faceRegion
                    );
                }
            }
        }
    }
    return nChanged;
}


// Mark all faces affected in any way by
// - removal of cells
// - removal of faces
// - removal of edges
// - removal of points
Foam::boolList Foam::removeFaces::getFacesAffected
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelList& facesToRemove,
    const labelHashSet& edgesToRemove,
    const labelHashSet& pointsToRemove
) const
{
    boolList affectedFace(mesh_.nFaces(), false);

    // Mark faces affected by removal of cells
    forAll(cellRegion, celli)
    {
        label region = cellRegion[celli];

        if (region != -1 && (celli != cellRegionMaster[region]))
        {
            const labelList& cFaces = mesh_.cells()[celli];

            forAll(cFaces, cFacei)
            {
                affectedFace[cFaces[cFacei]] = true;
            }
        }
    }

    // Mark faces affected by removal of face.
    forAll(facesToRemove, i)
    {
         affectedFace[facesToRemove[i]] = true;
    }

    //  Mark faces affected by removal of edges
    for (const label edgei : edgesToRemove)
    {
        const labelList& eFaces = mesh_.edgeFaces(edgei);

        forAll(eFaces, eFacei)
        {
            affectedFace[eFaces[eFacei]] = true;
        }
    }

    // Mark faces affected by removal of points
    for (const label pointi : pointsToRemove)
    {
        const labelList& pFaces = mesh_.pointFaces()[pointi];

        forAll(pFaces, pFacei)
        {
            affectedFace[pFaces[pFacei]] = true;
        }
    }
    return affectedFace;
}


void Foam::removeFaces::writeOBJ
(
    const indirectPrimitivePatch& fp,
    const fileName& fName
)
{
    OFstream str(fName);
    Pout<< "removeFaces::writeOBJ : Writing faces to file "
        << str.name() << endl;

    const pointField& localPoints = fp.localPoints();

    forAll(localPoints, i)
    {
        meshTools::writeOBJ(str, localPoints[i]);
    }

    const faceList& localFaces = fp.localFaces();

    forAll(localFaces, i)
    {
        const face& f = localFaces[i];

        str<< 'f';

        forAll(f, fp)
        {
            str<< ' ' << f[fp]+1;
        }
        str<< nl;
    }
}


// Get patch, zone info for facei
void Foam::removeFaces::getFaceInfo
(
    const label facei,

    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(facei))
    {
        patchID = mesh_.boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh_.faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}


// Return face with all pointsToRemove removed.
Foam::face Foam::removeFaces::filterFace
(
    const labelHashSet& pointsToRemove,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];

    labelList newFace(f.size(), -1);

    label newFp = 0;

    forAll(f, fp)
    {
        label vertI = f[fp];

        if (!pointsToRemove.found(vertI))
        {
            newFace[newFp++] = vertI;
        }
    }

    newFace.setSize(newFp);

    return face(newFace);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::removeFaces::removeFaces
(
    const polyMesh& mesh,
    const scalar minCos
)
:
    mesh_(mesh),
    minCos_(minCos)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Removing face connects cells. This function works out a consistent set of
// cell regions.
// - returns faces to remove. Can be extended with additional faces
//   (if owner would become neighbour)
// - sets cellRegion to -1 or to region number
// - regionMaster contains for every region number a master cell.
Foam::label Foam::removeFaces::compatibleRemoves
(
    const labelList& facesToRemove,
    labelList& cellRegion,
    labelList& regionMaster,
    labelList& newFacesToRemove
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    cellRegion.setSize(mesh_.nCells());
    cellRegion = -1;

    regionMaster.setSize(mesh_.nCells());
    regionMaster = -1;

    label nRegions = 0;

    forAll(facesToRemove, i)
    {
        label facei = facesToRemove[i];

        if (!mesh_.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Not internal face:" << facei << abort(FatalError);
        }


        label own = faceOwner[facei];
        label nei = faceNeighbour[facei];

        label region0 = cellRegion[own];
        label region1 = cellRegion[nei];

        if (region0 == -1)
        {
            if (region1 == -1)
            {
                // Create new region
                cellRegion[own] = nRegions;
                cellRegion[nei] = nRegions;

                // Make owner (lowest numbered!) the master of the region
                regionMaster[nRegions] = own;
                nRegions++;
            }
            else
            {
                // Add owner to neighbour region
                cellRegion[own] = region1;
                // See if owner becomes the master of the region
                regionMaster[region1] = min(own, regionMaster[region1]);
            }
        }
        else
        {
            if (region1 == -1)
            {
                // Add neighbour to owner region
                cellRegion[nei] = region0;
                // nei is higher numbered than own so guaranteed not lower
                // than master of region0.
            }
            else if (region0 != region1)
            {
                // Both have regions. Keep lowest numbered region and master.
                label freedRegion = -1;
                label keptRegion = -1;

                if (region0 < region1)
                {
                    changeCellRegion
                    (
                        nei,
                        region1,    // old region
                        region0,    // new region
                        cellRegion
                    );

                    keptRegion = region0;
                    freedRegion = region1;
                }
                else if (region1 < region0)
                {
                    changeCellRegion
                    (
                        own,
                        region0,    // old region
                        region1,    // new region
                        cellRegion
                    );

                    keptRegion = region1;
                    freedRegion = region0;
                }

                label master0 = regionMaster[region0];
                label master1 = regionMaster[region1];

                regionMaster[freedRegion] = -1;
                regionMaster[keptRegion] = min(master0, master1);
            }
        }
    }

    regionMaster.setSize(nRegions);


    // Various checks
    // - master is lowest numbered in any region
    // - regions have more than 1 cell
    {
        labelList nCells(regionMaster.size(), Zero);

        forAll(cellRegion, celli)
        {
            label r = cellRegion[celli];

            if (r != -1)
            {
                nCells[r]++;

                if (celli < regionMaster[r])
                {
                    FatalErrorInFunction
                        << "Not lowest numbered : cell:" << celli
                        << " region:" << r
                        << " regionmaster:" << regionMaster[r]
                        << abort(FatalError);
                }
            }
        }

        forAll(nCells, region)
        {
            if (nCells[region] == 1)
            {
                FatalErrorInFunction
                    << "Region " << region
                    << " has only " << nCells[region] << " cells in it"
                    << abort(FatalError);
            }
        }
    }


    // Count number of used regions
    label nUsedRegions = 0;

    forAll(regionMaster, i)
    {
        if (regionMaster[i] != -1)
        {
            nUsedRegions++;
        }
    }

    // Recreate facesToRemove to be consistent with the cellRegions.
    DynamicList<label> allFacesToRemove(facesToRemove.size());

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label own = faceOwner[facei];
        label nei = faceNeighbour[facei];

        if (cellRegion[own] != -1 && cellRegion[own] == cellRegion[nei])
        {
            // Both will become the same cell so add face to list of faces
            // to be removed.
            allFacesToRemove.append(facei);
        }
    }

    newFacesToRemove.transfer(allFacesToRemove);

    return nUsedRegions;
}


// ************************************************************************* //
