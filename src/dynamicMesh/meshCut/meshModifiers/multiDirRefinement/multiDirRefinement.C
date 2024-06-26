/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "multiDirRefinement.H"
#include "polyMesh.H"
#include "Time.H"
#include "undoableMeshCutter.H"
#include "hexCellLooper.H"
#include "geomCellLooper.H"
#include "directions.H"
#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "cellModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiDirRefinement, 0);
}


// * * * * * * * * * * * * * Private Statc Functions * * * * * * * * * * * * //

// Update refCells pattern for split cells. Note refCells is current
// list of cells to refine (these should all have been refined)
void Foam::multiDirRefinement::addCells
(
    const Map<label>& splitMap,
    List<refineCell>& refCells
)
{
    label newRefI = refCells.size();

    label oldSize = refCells.size();

    refCells.setSize(newRefI + splitMap.size());

    for (label refI = 0; refI < oldSize; refI++)
    {
        const refineCell& refCell = refCells[refI];

        const auto iter = splitMap.cfind(refCell.cellNo());

        if (!iter.good())
        {
            FatalErrorInFunction
                << "Problem : cannot find added cell for cell "
                << refCell.cellNo() << endl
                << abort(FatalError);
        }

        refCells[newRefI++] = refineCell(iter.val(), refCell.direction());
    }
}


// Update vectorField for all the new cells. Takes over value of
// original cell.
void Foam::multiDirRefinement::update
(
    const Map<label>& splitMap,
    vectorField& field
)
{
    field.setSize(field.size() + splitMap.size());

    forAllConstIters(splitMap, iter)
    {
        field[iter.val()] = field[iter.key()];
    }
}


// Add added cells to labelList
void Foam::multiDirRefinement::addCells
(
    const Map<label>& splitMap,
    labelList& labels
)
{
    label newCelli = labels.size();

    labels.setSize(labels.size() + splitMap.size());

    forAllConstIters(splitMap, iter)
    {
        labels[newCelli++] = iter.val();
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiDirRefinement::addCells
(
    const primitiveMesh& mesh,
    const Map<label>& splitMap
)
{
    // Construct inverse addressing: from new to original cell.
    labelList origCell(mesh.nCells(), -1);

    forAll(addedCells_, celli)
    {
        const labelList& added = addedCells_[celli];

        forAll(added, i)
        {
            label slave = added[i];

            if (origCell[slave] == -1)
            {
                origCell[slave] = celli;
            }
            else if (origCell[slave] != celli)
            {
                FatalErrorInFunction
                    << "Added cell " << slave << " has two different masters:"
                    << origCell[slave] << " , " << celli
                    << abort(FatalError);
            }
        }
    }


    forAllConstIters(splitMap, iter)
    {
        label masterI = iter.key();
        const label newCelli = iter.val();

        while (origCell[masterI] != -1 && origCell[masterI] != masterI)
        {
            masterI = origCell[masterI];
        }

        if (masterI >= addedCells_.size())
        {
            FatalErrorInFunction
                << "Map of added cells contains master cell " << masterI
                << " which is not a valid cell number" << endl
                << "This means that the mesh is not consistent with the"
                << " done refinement" << endl
                << "newCell:" << newCelli << abort(FatalError);
        }

        labelList& added = addedCells_[masterI];

        if (added.empty())
        {
            added.setSize(2);
            added[0] = masterI;
            added[1] = newCelli;
        }
        else if (!added.found(newCelli))
        {
            const label sz = added.size();
            added.setSize(sz + 1);
            added[sz] = newCelli;
        }
    }
}


Foam::labelList Foam::multiDirRefinement::splitOffHex(const primitiveMesh& mesh)
{
    const cellModel& hex = cellModel::ref(cellModel::HEX);

    const cellShapeList& cellShapes = mesh.cellShapes();

    // Split cellLabels_ into two lists.

    labelList nonHexLabels(cellLabels_.size());
    label nonHexI = 0;

    labelList hexLabels(cellLabels_.size());
    label hexI = 0;

    forAll(cellLabels_, i)
    {
        label celli = cellLabels_[i];

        if (cellShapes[celli].model() == hex)
        {
            hexLabels[hexI++] = celli;
        }
        else
        {
            nonHexLabels[nonHexI++] = celli;
        }
    }

    nonHexLabels.setSize(nonHexI);

    cellLabels_.transfer(nonHexLabels);

    hexLabels.setSize(hexI);

    return hexLabels;
}


void Foam::multiDirRefinement::refineHex8
(
    polyMesh& mesh,
    const labelList& hexCells,
    const bool writeMesh
)
{
    if (debug)
    {
        Pout<< "multiDirRefinement : Refining hexes " << hexCells.size()
            << endl;
    }

    hexRef8 hexRefiner
    (
        mesh,
        labelList(mesh.nCells(), Zero),    // cellLevel
        labelList(mesh.nPoints(), Zero),   // pointLevel
        refinementHistory
        (
            IOobject
            (
                "refinementHistory",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            List<refinementHistory::splitCell8>(),
            labelList(),
            false
        )                                   // refinement history
    );

    polyTopoChange meshMod(mesh);

    labelList consistentCells
    (
        hexRefiner.consistentRefinement
        (
            hexCells,
            true                  // buffer layer
        )
    );

    // Check that consistentRefinement hasn't added cells
    {
        // Create count 1 for original cells
        Map<label> hexCellSet(2*hexCells.size());
        for (const label celli : hexCells)
        {
            hexCellSet.insert(celli, 1);
        }

        // Increment count
        for (const label celli : consistentCells)
        {
            auto iter = hexCellSet.find(celli);

            if (iter.good())
            {
                iter.val() = 2;
            }
            else
            {
                FatalErrorInFunction
                    << "Resulting mesh would not satisfy 2:1 ratio"
                    << " when refining cell " << celli << abort(FatalError);
            }
        }

        // Check if all been visited (should always be since
        // consistentRefinement set up to extend set.
        forAllConstIters(hexCellSet, iter)
        {
            if (iter.val() != 2)
            {
                FatalErrorInFunction
                    << "Resulting mesh would not satisfy 2:1 ratio"
                    << " when refining cell " << iter.key()
                    << abort(FatalError);
            }
        }
    }


    hexRefiner.setRefinement(consistentCells, meshMod);

    // Change mesh, no inflation
    autoPtr<mapPolyMesh> morphMapPtr = meshMod.changeMesh(mesh, false, true);
    const mapPolyMesh& morphMap = morphMapPtr();

    if (morphMap.hasMotionPoints())
    {
        mesh.movePoints(morphMap.preMotionPoints());
    }

    if (writeMesh)
    {
        mesh.write();
    }

    if (debug)
    {
        Pout<< "multiDirRefinement : updated mesh at time "
            << mesh.time().timeName() << endl;
    }

    hexRefiner.updateMesh(morphMap);

    // Collect all cells originating from same old cell (original + 7 extra)

    forAll(consistentCells, i)
    {
        addedCells_[consistentCells[i]].setSize(8);
    }
    labelList nAddedCells(addedCells_.size(), Zero);

    const labelList& cellMap = morphMap.cellMap();

    forAll(cellMap, celli)
    {
        const label oldCelli = cellMap[celli];

        if (addedCells_[oldCelli].size())
        {
            addedCells_[oldCelli][nAddedCells[oldCelli]++] = celli;
        }
    }
}


void Foam::multiDirRefinement::refineAllDirs
(
    polyMesh& mesh,
    List<vectorField>& cellDirections,
    const cellLooper& cellWalker,
    undoableMeshCutter& cutter,
    const bool writeMesh
)
{
    // Iterator
    refinementIterator refiner(mesh, cutter, cellWalker, writeMesh);

    forAll(cellDirections, dirI)
    {
        if (debug)
        {
            Pout<< "multiDirRefinement : Refining " << cellLabels_.size()
                << " cells in direction " << dirI << endl
                << endl;
        }

        const vectorField& dirField = cellDirections[dirI];

        // Combine cell to be cut with direction to cut in.
        // If dirField is only one element use this for all cells.

        List<refineCell> refCells(cellLabels_.size());

        if (dirField.size() == 1)
        {
            // Uniform directions.
            if (debug)
            {
                Pout<< "multiDirRefinement : Uniform refinement:"
                    << dirField[0] << endl;
            }

            forAll(refCells, refI)
            {
                const label celli = cellLabels_[refI];

                refCells[refI] = refineCell(celli, dirField[0]);
            }
        }
        else
        {
            // Non uniform directions.
            forAll(refCells, refI)
            {
                const label celli = cellLabels_[refI];

                refCells[refI] = refineCell(celli, dirField[celli]);
            }
        }

        // Do refine mesh (multiple iterations). Remember added cells.
        Map<label> splitMap = refiner.setRefinement(refCells);

        // Update overall map of added cells
        addCells(mesh, splitMap);

        // Add added cells to list of cells to refine in next iteration
        addCells(splitMap, cellLabels_);

        // Update refinement direction for added cells.
        if (dirField.size() != 1)
        {
            forAll(cellDirections, i)
            {
                update(splitMap, cellDirections[i]);
            }
        }

        if (debug)
        {
            Pout<< "multiDirRefinement : Done refining direction " << dirI
                << " resulting in " << cellLabels_.size() << " cells" << nl
                << endl;
        }
    }
}


void Foam::multiDirRefinement::refineFromDict
(
    polyMesh& mesh,
    List<vectorField>& cellDirections,
    const dictionary& dict,
    const bool writeMesh
)
{
    // How to walk cell circumference.
    const bool pureGeomCut(dict.get<bool>("geometricCut"));

    autoPtr<cellLooper> cellWalker;
    if (pureGeomCut)
    {
        cellWalker.reset(new geomCellLooper(mesh));
    }
    else
    {
        cellWalker.reset(new hexCellLooper(mesh));
    }

    // Construct undoable refinement topology modifier.
    //Note: undoability switched off.
    // Might want to reconsider if needs to be possible. But then can always
    // use other constructor.
    undoableMeshCutter cutter(mesh, false);

    refineAllDirs(mesh, cellDirections, cellWalker(), cutter, writeMesh);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiDirRefinement::multiDirRefinement
(
    polyMesh& mesh,
    const labelList& cellLabels,        // list of cells to refine
    const dictionary& dict
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    const bool useHex(dict.get<bool>("useHexTopology"));

    const bool writeMesh(dict.get<bool>("writeMesh"));

    const wordList dirNames(dict.get<wordList>("directions"));

    if (useHex && dirNames.size() == 3)
    {
        // Filter out hexes from cellLabels_
        labelList hexCells(splitOffHex(mesh));

        refineHex8(mesh, hexCells, writeMesh);
    }

    if (returnReduceOr(cellLabels_.size()))
    {
        // Any cells to refine using meshCutter

        // Determine directions for every cell. Can either be uniform
        // (size = 1) or non-uniform (one vector per cell)
        directions cellDirections(mesh, dict);

        refineFromDict(mesh, cellDirections, dict, writeMesh);
    }
}


// Construct from dictionary and directions to refine.
Foam::multiDirRefinement::multiDirRefinement
(
    polyMesh& mesh,
    const labelList& cellLabels,        // list of cells to refine
    const List<vectorField>& cellDirs,  // Explicitly provided directions
    const dictionary& dict
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    const bool useHex(dict.get<bool>("useHexTopology"));

    const bool writeMesh(dict.get<bool>("writeMesh"));

    const wordList dirNames(dict.get<wordList>("directions"));

    if (useHex && dirNames.size() == 3)
    {
        // Filter out hexes from cellLabels_
        labelList hexCells(splitOffHex(mesh));

        refineHex8(mesh, hexCells, writeMesh);
    }

    if (returnReduceOr(cellLabels_.size()))
    {
        // Any cells to refine using meshCutter

        // Working copy of cells to refine
        List<vectorField> cellDirections(cellDirs);

        refineFromDict(mesh, cellDirections, dict, writeMesh);
    }
}


// Construct from components. Implies meshCutter since directions provided.
Foam::multiDirRefinement::multiDirRefinement
(
    polyMesh& mesh,
    undoableMeshCutter& cutter,     // actual mesh modifier
    const cellLooper& cellWalker,   // how to cut a single cell with
                                    // a plane
    const labelList& cellLabels,    // list of cells to refine
    const List<vectorField>& cellDirs,
    const bool writeMesh            // write intermediate meshes
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    // Working copy of cells to refine
    List<vectorField> cellDirections(cellDirs);

    refineAllDirs(mesh, cellDirections, cellWalker, cutter, writeMesh);
}


// ************************************************************************* //
