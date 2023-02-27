/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Alon Zameret, Noam Manaker Morag
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

#include "hpathRenumber.H"
#include "addToRunTimeSelectionTable.H"

#include "CircularBuffer.H"
#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hpathRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        hpathRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Local Details * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                 Class hpathRenumber::hpathFinder Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

// Class hpathFinder Declaration
class hpathRenumber::hpathFinder
{
    // Private Data

        // The input mesh
        const Foam::polyMesh& mesh;

        // Counter for the number of renumbered cells
        label nFoundCellCount;

        // Marks cells that have been added to the renumbering as 'true'
        Foam::bitSet bIsRenumbered;

        // For every data structure, I explain what it is used for and which
        // method is used to compute it

        // For every cell, a list of all its point-neighbouring cells - getMeshGraph()
        Foam::labelListList nMeshGraph;

        // For each cell its 'layer index': - getLayers()
        //      - cells in the same layer will have the same layer index
        Foam::labelList nCellLayerIndex;

        // For each cell its 'connected component index': - getConnectedComponents()
        //      - cells in the same connected component will have the same connected component index
        Foam::labelList nCellConnnectedComponentIndex;

        // For each cell its face-neighbours in the connected component - getConnnectedComponentGraph()
        Foam::labelListList nConnnectedComponentGraph;

        // Marks cells that have already been by BFS - getStartingCellInConnnectedComponent()
        Foam::bitSet bBFSFoundCell;

        // For each cell its point distance from the connected components start cell   - reorderDistFromStart()
        Foam::labelList nCellPointDistFromStart;
        Foam::labelList nCellFaceDistFromStart;

        // For each cell its DFS depth within the connected component     - findPath()
        Foam::labelList nDFSCellDepth;

        // For each cell its DFS parent within the connected component    - findPath()
        Foam::labelList nDFSParentCell;


    public: // Public methods

        //- Constructor
        explicit hpathFinder(const Foam::polyMesh& mesh);

        // Member Functions

        // Get Renumbering for the mesh
        void getRenumbering
        (
            Foam::labelList& cellOrder,
            bool bApplyLayerSeparation
        );

        // Returns the accuracy of the renumbering:
        // Accuracy is defined as the percentage of consecutive cells that are also face-neighbours in the mesh
        //  - Cells are face-neighbours if they have a common face
        float getAccuracy(const Foam::labelList& cellOrder) const;


    private: // Private methods

        // Initialize data structures for later use
        void initialize();

        // Given a cell and a face index, find its neighbour through the face
        // - If facing the boundary, returns -1
        // - Otherwise, returns FaceOwner/FaceNeighbour[nFaceIdx], the one that's different from nCellIdx
        label getNei(label nCellIdx, label nFaceIdx) const;

        // Creates the 'Mesh-Graph': for every cell, a list of cells that are point-neighbours with it in the mesh
        //  - Cells are point-neighbours if they have a common point
        void getMeshGraph();

        // Separates the mesh into layers: each cell has its layer saved in nCellLayerIndex
        // Also returns for every layer a list of all cells in it
        void getLayerSeparation(Foam::DynamicList<Foam::DynamicList<label>>& nCellsByLayer);

        // For every connected component, find the deepest cell and choose it as a starting cell
        // Returns a list with one starting cell per connected component
        void getStartingCells(Foam::DynamicList<label>& nStartingCells) const;

        // Once the starting cells have been found, this method does the actual layer separation
        void getLayers(const Foam::labelList& nStartingCellList, Foam::DynamicList<Foam::DynamicList<label>>& nCellsByLayer);

        // Renumber all cells in a layer
        void solveLayer(Foam::labelList nCellsInLayer, Foam::labelList& cellOrder);

        // Seperates cells into face-connected components
        // This is done using a general DFS algorithm
        void getConnectedComponents(const Foam::labelList& nCellList, Foam::DynamicList<Foam::DynamicList<label>>& nCellsByConnnectedComponent);

        // Find an approximate H-path through a connected component
        void solveConnectedComponent(const Foam::labelList& nCellsInConnectedComponent, Foam::labelList& cellOrder);

        // Finds for each cell in the connected component a list of its face-neighbours within the component
        void getConnnectedComponentGraph(const Foam::labelList& nCellsInConnectedComponent);

        // Finds a starting cell within the connected component
        label getStartingCellInConnnectedComponent(const Foam::labelList& nCellsInConnectedComponent);

        // This method reorders the cells in the given connected component based on their distance from nStartCell
        void reorderDistFromStart(label nStartCell, const Foam::labelList& nCellsInConnectedComponent);

        // Finds an H-path within the connected component and returns it in nResultHpath
        // H-path is guaranteed to start at nStartCell
        void findPath(label nStartCell, Foam::labelList& cellOrder);

        // Resets data structures for the cells that weren't found
        void resetCells(Foam::labelList& nCellList);
};

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hpathRenumber::hpathRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    applyLayerSeparation_
    (
        dict.optionalSubDict(typeName + "Coeffs")
            .getOrDefault("layered", true)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hpathRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    Info<< nl
        << "Starting Cell Renumbering: "
        << mesh.nCells() << " Cells, "
        << mesh.nFaces() << " Faces, "
        << mesh.nPoints() << " Points" << nl << nl;

    hpathRenumber::hpathFinder s(mesh);
    Foam::labelList cellOrder;
    s.getRenumbering(cellOrder, applyLayerSeparation_);

    std::cout << std::endl << "Found path with accuracy of "
        << std::fixed << std::setprecision(3)
        << s.getAccuracy(cellOrder) << "%" << std::endl;

    fprintf(stdout,"\n******************************************************\n\n\n");
    return cellOrder;
}


// * * * * * * * * * * * * hpathRenumber::hpathFinder  * * * * * * * * * * * //

// Public methods
Foam::hpathRenumber::hpathFinder::hpathFinder(const polyMesh& mesh)
:
    mesh(mesh)
{}


void Foam::hpathRenumber::hpathFinder::getRenumbering
(
    Foam::labelList& cellOrder,
    bool bApplyLayerSeparation
)
{
    // Find a renumbering for the entire mesh

    cellOrder.resize(mesh.nCells());
    // Counter for how many cells we have added to the renumbering so far
    nFoundCellCount = 0;

    // Initialize the data structures:
    Info<< "Initializing Data Structures" << endl;
    initialize();

    // Compute a graph to represent the mesh:
    //      - Cells in the mesh will be connected in the graph if they have a *common point*
    getMeshGraph();

    Foam::DynamicList<Foam::DynamicList<label>> nCellsByLayer;
    if (bApplyLayerSeparation)
    {
        getLayerSeparation(nCellsByLayer);
    }
    else
    {
        // If there is no layer separation, set the entire mesh as one 'layer'
        Foam::DynamicList<label> nAllCells(Foam::identity(mesh.nCells()));
        nCellsByLayer.append(nAllCells);
    }

    std::cout << "Beginning Hpath Computation" << std::endl;
    // Find H-path for each layer separately
    for (const Foam::labelList& nCellsInLayer : nCellsByLayer)
    {
        // solveLayer() will find a renumbering for all the cells in the layer
        // Path will be appended into cellOrder
        solveLayer(nCellsInLayer, cellOrder);
    }
}


float Foam::hpathRenumber::hpathFinder::getAccuracy
(
    const Foam::labelList& cellOrder
) const
{
    // Finding the number of "hits"
    // - A hit is a pair of consecutive cells in the renumbering that are also face-neighbours in the mesh
    // - Cells are face-neighbours if they have a common face
    label nHitCnt = 0;
    for (label nPathIdx = 0; nPathIdx < label(cellOrder.size()) - 1; nPathIdx++)
    {
        label nCurrCellIdx = cellOrder[nPathIdx];
        label nNextCellIdx = cellOrder[nPathIdx+1];

        // For every pair of consecutive cells, we search for a common face between them
        for (label nFaceIdx : mesh.cells()[nCurrCellIdx])
        {
            label nNeiIdx = getNei(nCurrCellIdx, nFaceIdx);
            // If there is a common face between them, we add a hit!
            if (nNeiIdx == nNextCellIdx)
            {
                nHitCnt++;
                break;
            }
        }
    }
    // The accuracy is the percentage of consecutive cells that were hits
    return 100.0f * float(nHitCnt) / float(label(cellOrder.size()) - 1);
}


// Private methods
void Foam::hpathRenumber::hpathFinder::initialize()
{
    // Data structure to keep track of cells already in the renumbering
    bIsRenumbered = Foam::bitSet(mesh.nCells(), false);

    // List used to seperate cells into layers
    // Saves for every cell its layer index
    nCellLayerIndex = Foam::labelList(mesh.nCells(), -1);

    // List used to seperate cells within the same layer into seperate connected components
    // Saves for every cell its connected component index
    nCellConnnectedComponentIndex = Foam::labelList(mesh.nCells(), -1);

    // Marks cells that have already been by BFS
    bBFSFoundCell = Foam::bitSet(mesh.nCells(), false);
    // Saves for every cell its face-neighbours within the same connected component
    nConnnectedComponentGraph = Foam::labelListList(mesh.nCells());

    // Saves for every cell within a layer its point-distance from the starting cell of that layer
    nCellPointDistFromStart = Foam::labelList(mesh.nCells(), -1);
    // Saves for every cell within a layer its face-distance from the starting cell of that layer
    nCellFaceDistFromStart = Foam::labelList(mesh.nCells(), -1);

    // For each cell its DFS depth within the connected component
    nDFSCellDepth = Foam::labelList(mesh.nCells(), -1);
    // For each cell its DFS parent within the connected component
    nDFSParentCell = Foam::labelList(mesh.nCells(), -1);

    // Now we can find the Mesh-Graph
    nMeshGraph = Foam::labelListList(mesh.nCells());
}


Foam::label Foam::hpathRenumber::hpathFinder::getNei(label nCell, label nFaceIdx) const
{
    // If the face has no neighbor, return -1
    if (nFaceIdx >= mesh.faceNeighbour().size())
    {
        return -1;
    }

    // Otherwise, the face connects 2 cells: the owner and the neighbor.
    // One of these should be nCell
    label nOwner = mesh.faceOwner()[nFaceIdx];
    label nNei = mesh.faceNeighbour()[nFaceIdx];

    // Find which of these two cells is the input cell, return the other one
    return (nOwner == nCell) ? nNei : nOwner;
}


void Foam::hpathRenumber::hpathFinder::getMeshGraph()
{
    // The purpose of the bitSet is to avoid adding the same cell as a neighbour more than once
    Foam::bitSet foundCells(mesh.nCells(), false);

    for (label nCellIdx = 0; nCellIdx < mesh.nCells(); nCellIdx++)
    {
        // Dynamic list of point-neighbours of the current cell
        DynamicList<label> nNeiList;

        for (label nPntIdx : mesh.cellPoints()[nCellIdx])
        {
            for (label nNeiCell : mesh.pointCells()[nPntIdx])
            {
                if (nNeiCell == nCellIdx) continue;
                if (foundCells[nNeiCell]) continue;

                // We have found a cell with a common point to the current cell
                // Therefore, we can now add it as a neighbour in the Mesh-Graph
                nNeiList.append(nNeiCell);

                foundCells[nNeiCell] = true;
            }
        }

        // Convert from DynamicList to labelList
        nMeshGraph[nCellIdx] = nNeiList;
        
        // Clear the bitSet before next iteration
        for (label neiCell : nMeshGraph[nCellIdx])
        {
            foundCells[neiCell] = false;
        }
    }
}


void Foam::hpathRenumber::hpathFinder::getLayerSeparation
(
    Foam::DynamicList<Foam::DynamicList<label>>& nCellsByLayer
)
{
    // We want to separate the mesh into layers:
    //      1) The mesh is separated into connected components                                          - getConnectedComponents()
    //      2) For each connected component we find the deepest cell and set it to be a starting cell   - getStartingCells()
    //      3) Each connected component is separated into layers using a BFS from the starting          - getLayers()

    std::cout << "Finding Layer Separation" << std::endl;

    // Step 1: Finding connected components
    Foam::labelList nAllCells(Foam::identity(mesh.nCells()));
    Foam::DynamicList<Foam::DynamicList<label>> nCellsByConnnectedComponent;
    getConnectedComponents(nAllCells, nCellsByConnnectedComponent);

    // Step 2: Finding starting cells
    // Finds for each connected component its deepest cell and returns them
    // The starting cells are returned through nStartingCellList
    // There will be exactly one starting cell per connected component
    Foam::DynamicList<label> nStartingCellList;
    getStartingCells(nStartingCellList);

    // Step 3: Layer separation
    // Layer separation is done in each connected component from the starting cell outwards
    getLayers(nStartingCellList, nCellsByLayer);

    label nLayerCnt = nCellsByLayer.size();
    std::cout << "Mesh Separated into " << nLayerCnt << " layers" << std::endl;

    // Reset connected component index list for solveLayer() to use
    nCellConnnectedComponentIndex = Foam::labelList(mesh.nCells(), -1);
}


void Foam::hpathRenumber::hpathFinder::getStartingCells
(
    Foam::DynamicList<label>& nStartingCells)
const
{
    // The starting cell in each connected component should be the 'deepest' cell in that connected component
    // A cells 'depth' is its point-distance from any boundary cell

    // The Algorithm:
    //      1) Find all boundary points
    //      2) Find all boundary cells
    //      3) Separate the boundary cells by connected component
    //      4) Find for each cell in the mesh its point-distance from the boundary
    //      5) Return through nStartingCells a list containing the deepest cell from each connected component

    // nCellDepthList will contain for each cell its depth within its ConnectedComponent
    Foam::labelList nCellDepthList(mesh.nCells(), -1);

    // Step 1: Finding all boundary points
    //      - A boundary point is a point on a boundary face
    Foam::bitSet bIsBoundaryPts(mesh.nPoints(), false);

    // Iterate over boundary faces and mark their points as boundary
    //  - Ignore boundary faces of type "empty"
    const polyBoundaryMesh& bndMesh = mesh.boundaryMesh();
    for (label nBndType = 0; nBndType < bndMesh.size(); ++nBndType)
    {
        // If boundary is of type "empty" - skip it
        if (bndMesh[nBndType].type().compare("empty") == 0) continue;

        labelRange range = bndMesh.patchRanges()[nBndType];
        // For every face in range set all of its points as boundary
        for (label nFaceIdx = range.min(); nFaceIdx <= range.max(); ++nFaceIdx)
        {
            for (label nPointIdx : mesh.faces()[nFaceIdx])
            {
                bIsBoundaryPts[nPointIdx] = true;
            }
        }
    }
    // Step 2: Finding all boundary cells
    //      - A boundary cell is a cell containing at least one boundary point
    Foam::bitSet bIsBoundaryCells(mesh.nCells(), false);

    // For every boundary point: mark all cells that have it as boundary cells
    for (label nPntIdx = 0; nPntIdx < mesh.nPoints(); nPntIdx++)
    {
        if (!bIsBoundaryPts[nPntIdx]) continue;
        // Here we make use of the nPntCellList we found when computing the Mesh-Graph
        for (label nCellIdx : mesh.pointCells()[nPntIdx])
        {
            bIsBoundaryCells[nCellIdx] = true;
        }
    }

    // Step 3: Separate the boundary cells based on which connected component they are a part of
    //      - Find for each connected component a list of all its boundary cells
    label nConnnectedComponentCnt = *std::max_element(nCellConnnectedComponentIndex.begin(), nCellConnnectedComponentIndex.end()) + 1;
    Foam::List<Foam::DynamicList<label>> nBoundaryCellsByConnnectedComponent(nConnnectedComponentCnt);
    for (label nCellIdx = 0; nCellIdx < mesh.nCells(); nCellIdx++)
    {
        if (bIsBoundaryCells[nCellIdx])
        {
            nBoundaryCellsByConnnectedComponent[nCellConnnectedComponentIndex[nCellIdx]].append(nCellIdx);
        }
    }

    // Steps 4+5:
    //      - Find each cells distance from the boundary by running a BFS algorithm from the boundary in each connected component
    //      - Along the way, find the deepest cell in each connected component and push them to the list

    for (Foam::labelList& nBndCells : nBoundaryCellsByConnnectedComponent)
    {
        // We want to find the maximum depth cell within the connected component
        label nMaxDepthCell = nBndCells[0];

        // Run a BFS algorithm from the boundary:
        //      - All boundary cells are pushed to the queue with depth 0
        Foam::CircularBuffer<label> nBfsQueue;
        for (label nCellIdx : nBndCells)
        {
            nCellDepthList[nCellIdx] = 0;
            nBfsQueue.push_back(nCellIdx);
        }

        while (!nBfsQueue.empty())
        {
            label nCurrCell = nBfsQueue.first();
            nBfsQueue.pop_front();

            // Check if current cell is the deeper than the maximum depth cell found so far. If it is, we update nMaxDepthCell
            if (nCellDepthList[nCurrCell] > nCellDepthList[nMaxDepthCell])
            {
                nMaxDepthCell = nCurrCell;
            }

            // The BFS algorithm needs to be based on point-neghbours, meaning by using the Mesh-Graph
            // Note that while cells in different connected components are never face-neighbours, but they may be point-neighbours (neighbours in the Mesh-Graph)
            // Therefore, when pushing all point-neighbouring cells we need to check that they are in the same connected component
            for (label nNeiCell : nMeshGraph[nCurrCell])
            {
                if (nCellDepthList[nNeiCell] != -1) continue;
                if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nCurrCell]) continue;

                nCellDepthList[nNeiCell] = nCellDepthList[nCurrCell] + 1;
                nBfsQueue.push_back(nNeiCell);
            }
        }

        // Finally, we can push the maximum depth cell we found
        nStartingCells.append(nMaxDepthCell);
    }
}


void Foam::hpathRenumber::hpathFinder::getLayers
(
    const Foam::labelList& nStartingCells,
    Foam::DynamicList<Foam::DynamicList<label>>& nCellsByLayer
)
{
    // Separates all cells in the mesh to layers
    // In each connected component:
    //      1) The starting cell is set as 'Layer 0'
    //      2) All cells that are point-neighbours of the starting cell are set as layer 1
    //      3) All cells that are point-neighbours of a cell in layer 1 are set as layer 2
    //      4) Repeat step (3) with increasing layer indices until all cells have been found
    // In this way, cells in each component are grouped in layers by their MINIMUM point-distance to their starting cell

    // Note: In reality, the same BFS is run for all components at once. Because components are connected, this is equivalent

    // Now we run BFS from starting cells
    Foam::CircularBuffer<label> nBfsQueue;

    // Push all starting cells to queue
    for (label nCellIdx : nStartingCells)
    {
        nCellLayerIndex[nCellIdx] = 0;
        nBfsQueue.push_back(nCellIdx);
    }

    // BFS algorithm will find for each cell its minimum distance in the Mesh-Graph from the corresponding starting cell
    while (!nBfsQueue.empty())
    {
        label nCurrCell = nBfsQueue.first();
        nBfsQueue.pop_front();

        label nLayer = nCellLayerIndex[nCurrCell];
        if (label(nCellsByLayer.size()) <= nCellLayerIndex[nCurrCell])
        {
            nCellsByLayer.append(Foam::labelList());
        }
        nCellsByLayer[nLayer].append(nCurrCell);

        for (label nNeiCell : nMeshGraph[nCurrCell])
        {
            // If neighbour is in a different connected component, ignore it
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nCurrCell]) continue;
            // If neighbour hasn't been visited yet, set it's layer and push it:
            if (nCellLayerIndex[nNeiCell] != -1) continue;
            nCellLayerIndex[nNeiCell] = nCellLayerIndex[nCurrCell] + 1;
            nBfsQueue.push_back(nNeiCell);
        }
    }
}


void Foam::hpathRenumber::hpathFinder::solveLayer
(
    Foam::labelList nCellsInLayer,
    Foam::labelList& cellOrder
)
{
    // Note: this method assumes all cells in the nCellsInLayer have the same 'layer index' in nCellLayerIndex

    // General rundown of the algorithm:
    //      1) Cells are separated into connected components
    //      2) Each connected component is solved separately and its path is added to the renumbering
    //      3) If there are still cells that have not been found, return to step (1)

    while (!nCellsInLayer.empty())
    {
        // Step 1: separate the cells into connected components
        //      - Connected components need to be connected by FACES (not points)
        Foam::DynamicList<Foam::DynamicList<label>> nCellsByConnnectedComponent;
        getConnectedComponents(nCellsInLayer, nCellsByConnnectedComponent);

        label nOrigFoundCellCount = nFoundCellCount;

        // Step 2: For each connected component we call getHpathinConnnectedComponent()
        for (const Foam::labelList& nCellsInConnectedComponent : nCellsByConnnectedComponent)
        {
            solveConnectedComponent(nCellsInConnectedComponent, cellOrder);
        }
        // Find how many cells were still not found
        label nRemainingCellCount = nCellsInLayer.size() - (nFoundCellCount - nOrigFoundCellCount);
        if (nRemainingCellCount == 0) break;

        // Step 3: Find all the cells in the layer that were not found yet
        Foam::labelList nRemainingCells(nRemainingCellCount);
        label nIndex = 0;
        for (label nCellIdx : nCellsInLayer)
        {
            if (!bIsRenumbered[nCellIdx])
            {
                nRemainingCells[nIndex++] = nCellIdx;
            }
        }

        resetCells(nRemainingCells);
        nCellsInLayer = nRemainingCells;
    }
}


void Foam::hpathRenumber::hpathFinder::getConnectedComponents
(
    const Foam::labelList& nCellList,
    Foam::DynamicList<Foam::DynamicList<label>>& nCellsByConnnectedComponent
)
{
    // This function separates cells in a certain layer into connected components
    // Each connected component will get a unique index
    //      - nCellsByConnnectedComponent[i] will contain a list of all cells with connected component index 'i'

    // We accomplish this using a generalized DFS algorithm:
    //      While there are still cells that haven't been found:
    //          1) Choose some cell in the layer that hasn't been found yet
    //          2) Find all cells connected to that cell (using a dfs stack)
    //          3) Set all of these cells as a new connected component

    // Notice: This method does not use the Mesh-Graph!
    //   - This is because the Mesh-Graph is POINT-connected, but we are searching for FACE-connected components

    // The index of the current connected component
    label nConnnectedComponentIndex = -1;

    for (label nCellIdx : nCellList)
    {
        if (nCellConnnectedComponentIndex[nCellIdx] != -1) continue;

        // This cell hasn't been found yet
        // So we set it to be the beginning of a new connected component
        // We increment nConnnectedComponentIndex, it is now the index of this new connected component
        nConnnectedComponentIndex++;
        nCellConnnectedComponentIndex[nCellIdx] = nConnnectedComponentIndex;
        nCellsByConnnectedComponent.append(Foam::labelList());

        Foam::CircularBuffer<label> nDfsStack;
        nDfsStack.push_back(nCellIdx);

        while(!nDfsStack.empty())
        {
            label nCurrCell = nDfsStack.last();
            nDfsStack.pop_back();

            nCellsByConnnectedComponent[nConnnectedComponentIndex].append(nCurrCell);

            // We push all face-neighbouring cells that:
            //      1) Are in the same layer
            //      2) Haven't been found yet
            for (label nFaceIdx : mesh.cells()[nCurrCell])
            {
                label nNeiCell = getNei(nCurrCell, nFaceIdx);
                if (nNeiCell < 0) continue;
                if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nCurrCell]) continue;
                if (nCellConnnectedComponentIndex[nNeiCell] != -1) continue;

                nCellConnnectedComponentIndex[nNeiCell] = nConnnectedComponentIndex;
                nDfsStack.push_back(nNeiCell);
            }
        }
    }
}


void Foam::hpathRenumber::hpathFinder::solveConnectedComponent
(
    const Foam::labelList& nCellsInConnectedComponent,
    Foam::labelList& cellOrder
)
{
    // For small cases, find path manually (1-2 cells)
    if (nCellsInConnectedComponent.size() <= 2) {
        for (label nCellIdx : nCellsInConnectedComponent) {
            cellOrder[nFoundCellCount++] = nCellIdx;
            bIsRenumbered[nCellIdx] = true;
        }
        return;
    }
    // Get a graph representing the connected component
    getConnnectedComponentGraph(nCellsInConnectedComponent);

    // Get the starting cell
    label nStartCell = getStartingCellInConnnectedComponent(nCellsInConnectedComponent);

    // Reorder each cells neighbours in the ConnnectedComponent-Graph in descending order by distance from start cell
    reorderDistFromStart(nStartCell, nCellsInConnectedComponent);

    // Get a path through the connected component, using the ConnnectedComponent-Graph
    findPath(nStartCell, cellOrder);
}


void Foam::hpathRenumber::hpathFinder::getConnnectedComponentGraph
(
    const Foam::labelList& nCellsInConnectedComponent
)
{
    // This method computes the ConnnectedComponent-Graph
    // Cells are connected in the ConnnectedComponent-Graph if:
    //      a) They are in the same layer
    //      b) They are face-connected in the mesh (different from Mesh-Graph - there it was point-connected)
    //          - note that this implies that they are also in the same connected component (hence the name)

    for (label nCellIdx : nCellsInConnectedComponent)
    {
        Foam::DynamicList<label> nNeiList;

        // For every face on the cell, find its neighbour through that face (if there is one)
        // If that neighbour:
        //      1) Hasn't been added to the reordering yet
        //      2) Is in the same layer
        //      3) Is in the same connected component
        // Then we add it as a neighbour in the ConnectedComponent-Graph
        for (label nFaceIdx : mesh.cells()[nCellIdx])
        {
            label nNeiCell = getNei(nCellIdx, nFaceIdx);
            if (nNeiCell < 0) continue;
            if (bIsRenumbered[nNeiCell]) continue;
            if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nCellIdx]) continue;

            nNeiList.append(nNeiCell);
        }

        nConnnectedComponentGraph[nCellIdx] = nNeiList;
    }
}


Foam::label Foam::hpathRenumber::hpathFinder::getStartingCellInConnnectedComponent
(
    const Foam::labelList& nCellsInConnectedComponent
)
{

    // Strategy for choosing a starting cell in the connected component:
    //      1) Start from an arbitrary cell within the connected component
    //      2) Find the furthest cell from it: choose it as the starting cell

    // The reasoning behind this strategy is as follows:
    //      - In many cases the cells given may be of the general shape of a long straight path
    //      - In these cases, by starting from the farthest cell we guarantee it will be at one of the edges of the path

    // It does not matter which cell we start from, so arbitrarly start from first cell in the connected component
    label nCellIdx = nCellsInConnectedComponent[0];

    // Find farthest cell from nCellIdx
    // This is done using a standard BFS algorithm
    Foam::CircularBuffer<label> nBfsQueue;

    bBFSFoundCell[nCellIdx] = true;
    nBfsQueue.push_back(nCellIdx);

    label nCurrCell = -1;

    while(!nBfsQueue.empty())
    {
        nCurrCell = nBfsQueue.first();
        nBfsQueue.pop_front();

        // Push all cells that have not yet been found by the BFS
        for (label nNeiCell : nConnnectedComponentGraph[nCurrCell])
        {
            if (bBFSFoundCell[nNeiCell]) continue;
            bBFSFoundCell[nNeiCell] = true;
            nBfsQueue.push_back(nNeiCell);
        }
    }

    // nCurrCell should now be the farthest cell from the starting cell in the connected component
    return nCurrCell;
}


void Foam::hpathRenumber::hpathFinder::reorderDistFromStart
(
    label nStartCell,
    const Foam::labelList& nCellsInConnectedComponent
)
{
    // findPath() tends to find much better results when each cell's neighbours are ordered in a specific way based on their face/point-distance to the starting cell

    // First, this method finds for each cell in the connected component its face-distance from the start cell
    // Secondly, this method finds for each cell in the connected component its point-distance from the start cell
    // Finally, this method reorders each cells neighbours in the ConnnectedComponent-Graph accordingly

    //  - Both of these are done using a BFS algorithm: the first using face-neigbours and the second using point-neighbours

    Foam::CircularBuffer<label> nBfsQueue;

    // Find the face-distance of all cells in the connected component from the starting cell
    nCellFaceDistFromStart[nStartCell] = 0;
    // Push starting cell to queue
    nBfsQueue.push_back(nStartCell);

    while (!nBfsQueue.empty())
    {
        label nCurrCell = nBfsQueue.first();
        nBfsQueue.pop_front();

        // Push all face-neighbours that haven't already had their face distance found
        // Face-neighbours are saved in the ConnnectedComponent-Graph
        for (label nNeiCell : nConnnectedComponentGraph[nCurrCell])
        {
            if (nCellFaceDistFromStart[nNeiCell] != -1) continue;

            nCellFaceDistFromStart[nNeiCell] = nCellFaceDistFromStart[nCurrCell] + 1;
            nBfsQueue.push_back(nNeiCell);
        }
    }
    // Find the *point*-distance of all cells in the connected component from the starting cell
    nCellPointDistFromStart[nStartCell] = 0;
    // Push starting cell to queue
    nBfsQueue.push_back(nStartCell);

    while (!nBfsQueue.empty())
    {
        label nCurrCell = nBfsQueue.first();
        nBfsQueue.pop_front();

        // Push all the point-neighbour that:
        //      1) Haven't already been pushed by the BFS previously
        //      2) Are in the same layer as the Starting Cell
        //      3) Are in the same connected component as the Starting Cell
        // Point neighbours are saved in the Mesh-Graph
        for (label nNeiCell : nMeshGraph[nCurrCell])
        {
            if (bIsRenumbered[nNeiCell]) continue;
            if (nCellPointDistFromStart[nNeiCell] != -1) continue;
            if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nStartCell]) continue;
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nStartCell]) continue;

            nCellPointDistFromStart[nNeiCell] = nCellPointDistFromStart[nCurrCell] + 1;
            nBfsQueue.push_back(nNeiCell);
        }
    }

    for (label nCellIdx : nCellsInConnectedComponent)
    {
        // Sorting of neighbours is done in two levels:
        //      1) Neigbours are sorted *descending*-order based on their *point*-distance from the starting cell
        //      2) Neigbours with the same *point*-distance are sorted in *ascending*-order based on their *face*-distance from the starting cell
        std::sort
        (
            nConnnectedComponentGraph[nCellIdx].begin(),
            nConnnectedComponentGraph[nCellIdx].end(),
            [this](label i, label j) {
            if (nCellPointDistFromStart[i] != nCellPointDistFromStart[j])
                return nCellPointDistFromStart[i] > nCellPointDistFromStart[j];
            else
                return nCellFaceDistFromStart[i] < nCellFaceDistFromStart[j];
            }
        );
    }
}


void Foam::hpathRenumber::hpathFinder::findPath
(
    label nStartCell, Foam::labelList& cellOrder
)
{
    // Tries to find the furthest cell from the starting cell within the ConnectedComponent-Graph
    // This is done using a DFS algorithm:
    //      - Run DFS from starting cell
    //      - Return path to deepest cell found by DFS

    // Run a standard DFS algorithm from starting cell
    Foam::CircularBuffer<label> nDfsStack;
    nDfsStack.push_back(nStartCell);
    nDFSParentCell[nStartCell] = nStartCell;

    // Used for finding and returning the best path
    label nBestCell = nStartCell;

    while (!nDfsStack.empty())
    {
        label nCurrCell = nDfsStack.last();
        nDfsStack.pop_back();

        // If the cell has already been popped previously, we can skip it
        if (nDFSCellDepth[nCurrCell] != -1) continue;

        // Otherwise, we update its depth value. This means the cell will never be PUSHED again
        nDFSCellDepth[nCurrCell] = nDFSCellDepth[nDFSParentCell[nCurrCell]] + 1;

        // If this is the new deepest cell, update best
        if (nDFSCellDepth[nCurrCell] > nDFSCellDepth[nBestCell]) {
            nBestCell = nCurrCell;
        }

        // Cells will be popped from the stack in reverse order from how we pushed them
        // By iterating over the neighbors in reverse order, cells will be popped in the correct order
        for (label nNeiIdx = nConnnectedComponentGraph[nCurrCell].size() - 1; nNeiIdx >= 0; nNeiIdx--)
        {
            label nNeiCell = nConnnectedComponentGraph[nCurrCell][nNeiIdx];
            if (nDFSCellDepth[nNeiCell] == -1)
            {
                // Notice we only update the nCellDepth value of a cell when we POP it from the stack
                // This means some cells may be pushed many times, but they will only be popped once
                nDFSParentCell[nNeiCell] = nCurrCell;
                nDfsStack.push_back(nNeiCell);
            }
        }
    }

    label nCellIdx = nBestCell;
    label nPrevCell = -1;
    while (nCellIdx != nPrevCell)
    {
        cellOrder[nFoundCellCount++] = nCellIdx;
        bIsRenumbered[nCellIdx] = true;

        // The first cell in the path is always its own parent
        nPrevCell = nCellIdx;
        nCellIdx = nDFSParentCell[nCellIdx];
    }
}


void Foam::hpathRenumber::hpathFinder::resetCells
(
    Foam::labelList& nCellList
)
{
    // Reset the data structures for the cells that were not found in the renumbering
    for (label nCellIdx : nCellList)
    {
        nCellConnnectedComponentIndex[nCellIdx] = -1;
        nCellPointDistFromStart[nCellIdx] = -1;
        nCellFaceDistFromStart[nCellIdx] = -1;
        bBFSFoundCell[nCellIdx] = false;
        nDFSCellDepth[nCellIdx] = -1;
        nDFSParentCell[nCellIdx] = -1;
    }
}


// ************************************************************************* //
