/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "removeCells.H"
#include "bitSet.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "polyRemoveCell.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(removeCells, 0);
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Increase count (usage) of elements of list
inline void incrCount(const Foam::labelUList& list, Foam::labelList& counter)
{
    for (auto idx : list)
    {
        ++counter[idx];
    }
}

// Decrease count (usage) of elements of list
inline void decrCount(const Foam::labelUList& list, Foam::labelList& counter)
{
    for (auto idx : list)
    {
        --counter[idx];
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::removeCells::removeCells(const polyMesh& mesh, const bool syncPar)
:
    mesh_(mesh),
    syncPar_(syncPar)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::removeCells::getExposedFaces
(
    const bitSet& removedCell
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Count cells using face.
    labelList nCellsUsingFace(mesh_.nFaces(), Zero);

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        const label own = faceOwner[facei];
        const label nei = faceNeighbour[facei];

        if (!removedCell[own])
        {
            ++nCellsUsingFace[facei];
        }
        if (!removedCell[nei])
        {
            ++nCellsUsingFace[facei];
        }
    }

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        const label own = faceOwner[facei];

        if (!removedCell[own])
        {
            ++nCellsUsingFace[facei];
        }
    }

    // Coupled faces: add number of cells using face across couple.
    {
        // Note cyclics done always, parallel bits only done if syncPar_

        SubList<label> bndValues
        (
            nCellsUsingFace,
            mesh_.nBoundaryFaces(),
            mesh_.nInternalFaces()
        );

        syncTools::syncBoundaryFaceList
        (
            mesh_,
            bndValues,
            plusEqOp<label>(),
            mapDistribute::transform(),
            syncPar_
        );
    }

    // Now nCellsUsingFace:
    // 0 : internal face whose both cells get deleted
    //     boundary face whose all cells get deleted
    // 1 : internal face that gets exposed
    //     unaffected (uncoupled) boundary face
    //     coupled boundary face that gets exposed ('uncoupled')
    // 2 : unaffected internal face
    //     unaffected coupled boundary face

    DynamicList<label> exposedFaces(mesh_.nFaces()/10);

    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (nCellsUsingFace[facei] == 1)
        {
            exposedFaces.append(facei);
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                const label own = faceOwner[facei];

                if (nCellsUsingFace[facei] == 1 && !removedCell[own])
                {
                    // My owner not removed but other side is so has to become
                    // normal, uncoupled, boundary face
                    exposedFaces.append(facei);
                }

                ++facei;
            }
        }
    }

    return labelList(std::move(exposedFaces));
}


void Foam::removeCells::setRefinement
(
    const bitSet& removedCell,
    const labelUList& exposedFaceLabels,
    const labelUList& exposedPatchIDs,
    polyTopoChange& meshMod
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (exposedFaceLabels.size() != exposedPatchIDs.size())
    {
        FatalErrorInFunction
            << "Size of exposedFaceLabels " << exposedFaceLabels.size()
            << " differs from size of exposedPatchIDs "
            << exposedPatchIDs.size()
            << abort(FatalError);
    }

    // List of new patchIDs
    labelList newPatchID(mesh_.nFaces(), -1);

    forAll(exposedFaceLabels, i)
    {
        const label facei = exposedFaceLabels[i];
        const label patchi = exposedPatchIDs[i];

        if (patchi < 0 || patchi >= patches.size())
        {
            FatalErrorInFunction
                << "Invalid patch " << patchi
                << " for exposed face " << facei << nl
                << "Valid patches 0.." << patches.size()-1
                << abort(FatalError);
        }

        if (patches[patchi].coupled())
        {
            FatalErrorInFunction
                << "Trying to put exposed face " << facei
                << " into a coupled patch : " << patches[patchi].name()
                << nl
                << "This is illegal."
                << abort(FatalError);
        }

        newPatchID[facei] = patchi;
    }


    // Walk all the cells mentioned for removal
    for (const label celli : removedCell)
    {
        //Pout<< "Removing cell " << celli
        //    << " cc:" << mesh_.cellCentres()[celli] << endl;

        meshMod.setAction(polyRemoveCell(celli));
    }


    // Remove faces that are no longer used. Modify faces that
    // are used by one cell only.

    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const faceZoneMesh& faceZones = mesh_.faceZones();

    // Count starting number of faces using each point.
    // Update whenever removing a face.
    labelList nFacesUsingPoint(mesh_.nPoints(), Zero);

    for (const face& f : faces)
    {
        incrCount(f, nFacesUsingPoint);
    }


    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        const face& f = faces[facei];
        const label own = faceOwner[facei];
        const label nei = faceNeighbour[facei];

        if (removedCell[own])
        {
            if (removedCell[nei])
            {
                // Face no longer used
                //Pout<< "Removing internal face " << facei
                //    << " fc:" << mesh_.faceCentres()[facei] << endl;

                meshMod.setAction(polyRemoveFace(facei));
                decrCount(f, nFacesUsingPoint);
            }
            else
            {
                if (newPatchID[facei] == -1)
                {
                    FatalErrorInFunction
                        << "No patchID provided for exposed face " << facei
                        << " on cell " << nei << nl
                        << "Did you provide patch IDs for all exposed faces?"
                        << abort(FatalError);
                }

                // nei is remaining cell. Facei becomes external cell

                const label zoneID = faceZones.whichZone(facei);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    // Note: we reverse the owner/neighbour of the face
                    // so should also select the other side of the zone
                    zoneFlip = !fZone.flipMap()[fZone.whichFace(facei)];
                }

                //Pout<< "Putting exposed internal face " << facei
                //    << " fc:" << mesh_.faceCentres()[facei]
                //    << " into patch " << newPatchID[facei] << endl;

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        f.reverseFace(),        // modified face
                        facei,                  // label of face being modified
                        nei,                    // owner
                        -1,                     // neighbour
                        true,                   // face flip
                        newPatchID[facei],      // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
        else if (removedCell[nei])
        {
            if (newPatchID[facei] == -1)
            {
                FatalErrorInFunction
                    << "No patchID provided for exposed face " << facei
                    << " on cell " << own << nl
                    << "Did you provide patch IDs for all exposed faces?"
                    << abort(FatalError);
            }

            //Pout<< "Putting exposed internal face " << facei
            //    << " fc:" << mesh_.faceCentres()[facei]
            //    << " into patch " << newPatchID[facei] << endl;

            // own is remaining cell. Facei becomes external cell.
            const label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                      // modified face
                    facei,                  // label of face being modified
                    own,                    // owner
                    -1,                     // neighbour
                    false,                  // face flip
                    newPatchID[facei],      // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if (newPatchID[facei] != -1)
                {
                    //Pout<< "Putting uncoupled coupled face " << facei
                    //    << " fc:" << mesh_.faceCentres()[facei]
                    //    << " into patch " << newPatchID[facei] << endl;

                    const label zoneID = faceZones.whichZone(facei);
                    bool zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            faces[facei],           // modified face
                            facei,                  // label of face
                            faceOwner[facei],       // owner
                            -1,                     // neighbour
                            false,                  // face flip
                            newPatchID[facei],      // patch for face
                            false,                  // remove from zone
                            zoneID,                 // zone for face
                            zoneFlip                // face flip in zone
                        )
                    );
                }
                else if (removedCell[faceOwner[facei]])
                {
                    // Face no longer used
                    //Pout<< "Removing boundary face " << facei
                    //    << " fc:" << mesh_.faceCentres()[facei]
                    //    << endl;

                    meshMod.setAction(polyRemoveFace(facei));
                    decrCount(faces[facei], nFacesUsingPoint);
                }

                ++facei;
            }
        }
        else
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if (newPatchID[facei] != -1)
                {
                    FatalErrorInFunction
                        << "new patchID provided for boundary face " << facei
                        << " even though it is not on a coupled face."
                        << abort(FatalError);
                }

                if (removedCell[faceOwner[facei]])
                {
                    // Face no longer used
                    //Pout<< "Removing boundary face " << facei
                    //    << " fc:" << mesh_.faceCentres()[facei]
                    //    << endl;

                    meshMod.setAction(polyRemoveFace(facei));
                    decrCount(faces[facei], nFacesUsingPoint);
                }

                ++facei;
            }
        }
    }


    // Remove points that are no longer used.
    // Loop rewritten to not use pointFaces.

    forAll(nFacesUsingPoint, pointi)
    {
        if (nFacesUsingPoint[pointi] == 0)
        {
            //Pout<< "Removing unused point " << pointi
            //    << " at:" << mesh_.points()[pointi] << endl;

            meshMod.setAction(polyRemovePoint(pointi));
        }
        else if (nFacesUsingPoint[pointi] == 1)
        {
            WarningInFunction
                << "point " << pointi << " at coordinate "
                << mesh_.points()[pointi]
                << " is only used by 1 face after removing cells."
                << " This probably results in an illegal mesh."
                << endl;
        }
    }
}


Foam::labelList Foam::removeCells::getExposedFaces
(
    const labelUList& cellsToRemove
) const
{
    bitSet removeCell(mesh_.nCells(), cellsToRemove);

    return getExposedFaces(removeCell);
}


void Foam::removeCells::setRefinement
(
    const labelUList& cellsToRemove,
    const labelUList& exposedFaceLabels,
    const labelUList& exposedPatchIDs,
    polyTopoChange& meshMod
) const
{
    bitSet removedCell(mesh_.nCells(), cellsToRemove);

    setRefinement
    (
        removedCell,
        exposedFaceLabels,
        exposedPatchIDs,
        meshMod
    );
}


// ************************************************************************* //
