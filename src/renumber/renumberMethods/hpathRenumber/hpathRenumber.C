/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include <queue>
#include <stack>
#include <iomanip>
#include <numeric>

#include "hpathRenumber.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hpathRenumber::hpathRenumber(const dictionary& renumberDict) :
    renumberMethod(renumberDict),
    m_bApplyLayerSeparation(renumberDict.optionalSubDict(typeName + "Coeffs").getOrDefault("layered", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::labelList Foam::hpathRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    fprintf(stdout,"\n\n\n******************************************************\n\n");
    fprintf(stdout, "Starting Cell Renumbering: %d Cells, %d Faces, %d Points\n\n", mesh.nCells(), mesh.nFaces(), mesh.nPoints());

    hpathFinder s(mesh);
    Foam::labelList cellOrder;
    s.getRenumbering(cellOrder, m_bApplyLayerSeparation);

    std::cout << std::endl << "Found path with accuracy of " << std::fixed << std::setprecision(3) << s.getAccuracy(cellOrder) << "%" << std::endl;

    fprintf(stdout,"\n******************************************************\n\n\n");
    return cellOrder;
}


// * * * * * * * * * * * * * * * hpathRenumber::dynamicMarker  * * * * * * * * * * * * * * * //

// This class supports marking and unmarking small sets of indices quickly
// Used for efficient construction of Mesh Graph in getMeshGraph()

// At initialization all items begin unmarked (false)
Foam::hpathRenumber::dynamicMarker::dynamicMarker(int n) {
    bMarker.assign(n, false);
}


// When we mark an index, we also push it to the stack
bool Foam::hpathRenumber::dynamicMarker::mark(int i) {
    if (bMarker[i]) return false;
    bMarker[i] = true;
    nMarkedStack.push(i);
    return true;
}

// When we want to clear all marks, we unmark every item in the stack
// This takes amortized-time O(1), because we can only pop as many items as we have pushed
void Foam::hpathRenumber::dynamicMarker::clear() {
    while(!nMarkedStack.empty()) {
        bMarker[nMarkedStack.top()] = false;
        nMarkedStack.pop();
    }
}


// * * * * * * * * * * * * * * * hpathRenumber::hpathFinder * * * * * * * * * * * * * * * //

// Public methods

Foam::hpathRenumber::hpathFinder::hpathFinder(const polyMesh& mesh) : mesh(mesh), nCellCount(mesh.nCells()) {}


void Foam::hpathRenumber::hpathFinder::getRenumbering(Foam::labelList& cellOrder, bool bApplyLayerSeparation) {

    // Find a renumbering for the entire mesh

    cellOrder.resize(nCellCount);
    // Counter for how many cells we have added to the renumbering so far
    nFoundCellCount = 0;

    // Initialize the data structures:
    std::cout << "Initializing Data Structures" << std::endl;
    initialize();
    
    // Compute a graph to represent the mesh:
    //      - Cells in the mesh will be connected in the graph if they have a *common point*
    getMeshGraph();

    std::vector<std::vector<int>> nCellsByLayer;
    if (bApplyLayerSeparation) {
        getLayerSeparation(nCellsByLayer);
    }
    else {
        // If there is no layer separation, set the entire mesh as one 'layer'
        nCellsByLayer.emplace_back(nCellCount,-1);
        std::iota(nCellsByLayer[0].begin(), nCellsByLayer[0].end(), 0);
    }

    std::cout << "Beginning Hpath Computation" << std::endl;
    // Find H-path for each layer separately
    for (const std::vector<int>& nCellsInLayer : nCellsByLayer)
    {
        // solveLayer() will find a renumbering for all the cells in the layer
        // Path will be appended into cellOrder
        solveLayer(nCellsInLayer, cellOrder);
    }
}


float Foam::hpathRenumber::hpathFinder::getAccuracy(const Foam::labelList& cellOrder) const {
    // Finding the number of "hits"
    // - A hit is a pair of consecutive cells in the renumbering that are also face-neighbours in the mesh
    // - Cells are face-neighbours if they have a common face
    int nHitCnt = 0;
    for (int nPathIdx = 0; nPathIdx < int(cellOrder.size()) - 1; nPathIdx++)
    {
        int nCurrCellIdx = cellOrder[nPathIdx];
        int nNextCellIdx = cellOrder[nPathIdx+1];

        // For every pair of consecutive cells, we search for a common face between them
        for (int nFaceIdx : mesh.cells()[nCurrCellIdx])
        {
            int nNeiIdx = getNei(nCurrCellIdx, nFaceIdx);
            // If there is a common face between them, we add a hit!
            if (nNeiIdx == nNextCellIdx) {
                nHitCnt++;
                break;
            }
        }
    }
    // The accuracy is the percentage of consecutive cells that were hits
    return 100.0f * float(nHitCnt) / float(int(cellOrder.size()) - 1);
}

// Private methods

void Foam::hpathRenumber::hpathFinder::initialize() {
    
    // Data structure to keep track of cells already in the renumbering
    bIsRenumbered.assign(nCellCount, false);

    // List used to seperate cells int layers
    // Saves for every cell its layer index
    nCellLayerIndex.assign(nCellCount, -1);

    // List used to seperate cells within the same layer into seperate connected components
    // Saves for every cell its connected component index
    nCellConnnectedComponentIndex.assign(nCellCount, -1);

    // Marks cells that have already been by BFS
    bBFSFoundCell.assign(nCellCount, false);
    // Saves for every cell its face-neighbours within the same connected component
    nConnnectedComponentGraph.assign(nCellCount, std::vector<int>());

    // Saves for every cell within a layer its point-distance from the starting cell of that layer
    nCellPointDistFromStart.assign(nCellCount, -1);
    // Saves for every cell within a layer its face-distance from the starting cell of that layer
    nCellFaceDistFromStart.assign(nCellCount, -1);

    // For each cell its DFS depth within the connected component
    nDFSCellDepth.assign(nCellCount, -1);
    // For each cell its DFS parent within the connected component
    nDFSParentCell.assign(nCellCount, -1);

    // Now we can find the Mesh-Graph
    nMeshGraph.assign(nCellCount, std::vector<int>());
}


int Foam::hpathRenumber::hpathFinder::getNei(int nCell, int nFaceIdx) const {
    // If the face has no neighbor, return -1
    if(nFaceIdx >= mesh.faceNeighbour().size()) 
        return -1;
    
    // Otherwise, the face connects 2 cells: the owner and the neighbor. One of these should be nCell
    int nOwner = mesh.faceOwner()[nFaceIdx];
    int nNei = mesh.faceNeighbour()[nFaceIdx];

    // Find which of these two cells is the input cell, return the other one
    return (nOwner == nCell) ? nNei : nOwner;
}

void Foam::hpathRenumber::hpathFinder::getMeshGraph() {

    // First, we need to compute:
    //      1) For every point a list of cells it is part of  - nPntCellList
    //      2) For every cell a list of points on it          - nCellPntList
    getPointCellLists();

    // The purpose of the dynamic marker is to avoid checking the same cell many times in one iteration
    dynamicMarker bCellMarker(nCellCount);
    for (int nCellIdx = 0; nCellIdx < nCellCount; nCellIdx++) {
        for (int nPntIdx : nCellPntList[nCellIdx]) {
            for (int nNeiCell : nPntCellList[nPntIdx]) {
                if (nNeiCell == nCellIdx) continue;
                if (!bCellMarker.mark(nNeiCell)) continue;

                // We have found a cell with a common point to the current cell
                // Therefore, we can now add it as a neighbour in the Mesh-Graph
                nMeshGraph[nCellIdx].push_back(nNeiCell);
            }
        }
        bCellMarker.clear();
    }
}


void Foam::hpathRenumber::hpathFinder::getPointCellLists() {

    // Find:
    //  - For each cell in the mesh a list of all points on it
    nCellPntList.assign(mesh.nCells(),  std::vector<int>());
    //  - For every point in the mesh a list of all cells it is part of
    nPntCellList.assign(mesh.nPoints(), std::vector<int>());

    // The purpose of the dynamic marker is to avoid pushing the same point many times in one iteration
    dynamicMarker bPointMarker(mesh.nPoints());
    for (int nCellIdx = 0; nCellIdx < mesh.nCells(); nCellIdx++)
    {
        for(int nFaceIdx : mesh.cells()[nCellIdx]) {
            for(int nPntIdx : mesh.faces()[nFaceIdx]) {
                // If we already found nPntIdx for this cell, continue
                if (!bPointMarker.mark(nPntIdx)) continue;
                
                // nPntIdx is on nCellIdx, so we push them to each-others lists
                nCellPntList[nCellIdx].push_back(nPntIdx);
                nPntCellList[nPntIdx].push_back(nCellIdx);
            }
        }
        // Clear the marked points for the next cell
        bPointMarker.clear();
    }
}

void Foam::hpathRenumber::hpathFinder::getLayerSeparation(std::vector<std::vector<int>>& nCellsByLayer)
{
    // We want to separate the mesh into layers:
    //      1) The mesh is separated into connected components                                          - getConnectedComponents()
    //      2) For each connected component we find the deepest cell and set it to be a starting cell   - getStartingCells()
    //      3) Each connected component is separated into layers using a BFS from the starting          - getLayers()

    std::cout << "Finding Layer Separation" << std::endl;

    // Step 1: Finding connected components
    std::vector<int> nAllCells(nCellCount);
    std::iota(nAllCells.begin(), nAllCells.end(), 0);
    
    std::vector<std::vector<int>> nCellsByConnnectedComponent;
    getConnectedComponents(nAllCells, nCellsByConnnectedComponent);

    // Step 2: Finding starting cells
    // Finds for each connected component its deepest cell and returns them
    // The starting cells are returned through nStartingCellList
    // There will be exactly one starting cell per connected component
    std::vector<int> nStartingCellList;
    getStartingCells(nStartingCellList);

    // Step 3: Layer separation
    // Layer separation is done in each connected component from the starting cell outwards
    getLayers(nStartingCellList, nCellsByLayer);
    
    int nLayerCnt = nCellsByLayer.size();
    std::cout << "Mesh Separated into " << nLayerCnt << " layers" << std::endl;
    
    // Reset connected component index list for solveLayer() to use
    nCellConnnectedComponentIndex.assign(nCellCount, -1);
}

void Foam::hpathRenumber::hpathFinder::getStartingCells(std::vector<int>& nStartingCells) const {

    // The starting cell in each connected component should be the 'deepest' cell in that connected component
    // A cells 'depth' is its point-distance from any boundary cell
    
    // The Algorithm:
    //      1) Find all boundary points
    //      2) Find all boundary cells
    //      3) Separate the boundary cells by connected component
    //      4) Find for each cell in the mesh its point-distance from the boundary
    //      5) Return through nStartingCells a list containing the deepest cell from each connected component

    // nCellDepthList will contain for each cell its depth within its ConnectedComponent
    std::vector<int> nCellDepthList(nCellCount, -1);

    // Step 1: Finding all boundary points
    //      - A boundary point is a point on a boundary face
    std::vector<int> bIsBoundaryPts(mesh.nPoints(), false);

    // Iterate over boundary faces and mark their points as boundary
    //  - Ignore boundary faces of type "empty"
    const polyBoundaryMesh& bndMesh = mesh.boundaryMesh();
    for (int nBndType = 0; nBndType < bndMesh.size(); ++nBndType)
    {
        // If boundary is of type "empty" - skip it
        if (bndMesh[nBndType].type().compare("empty") == 0) continue;

        labelRange range = bndMesh.patchRanges()[nBndType];
        // For every face in range set all of its points as boundary
        for (int nFaceIdx = range.min(); nFaceIdx <= range.max(); ++nFaceIdx) {
            for(int nPointIdx : mesh.faces()[nFaceIdx]) {
                bIsBoundaryPts[nPointIdx] = true;
            }
        }
    }
    
    // Step 2: Finding all boundary cells
    //      - A boundary cell is a cell containing at least one boundary point
    std::vector<bool> bIsBoundaryCells(nCellCount, false);

    // For every boundary point: mark all cells that have it as boundary cells
    for (int nPntIdx = 0; nPntIdx < mesh.nPoints(); nPntIdx++)
    {
        if (!bIsBoundaryPts[nPntIdx]) continue;
        // Here we make use of the nPntCellList we found when computing the Mesh-Graph
        for (int nCellIdx : nPntCellList[nPntIdx])
            bIsBoundaryCells[nCellIdx] = true;
    }

    // Step 3: Separate the boundary cells based on which connected component they are a part of
    //      - Find for each connected component a list of all its boundary cells
    int nConnnectedComponentCnt = *std::max_element(nCellConnnectedComponentIndex.begin(), nCellConnnectedComponentIndex.end()) + 1;
    std::vector<std::vector<int>> nBoundaryCellsByConnnectedComponent(nConnnectedComponentCnt, std::vector<int>());
    for (int nCellIdx = 0; nCellIdx < nCellCount; nCellIdx++)
    {
        if (bIsBoundaryCells[nCellIdx]) {
            nBoundaryCellsByConnnectedComponent[nCellConnnectedComponentIndex[nCellIdx]].push_back(nCellIdx);
        }
    }

    // Steps 4+5:
    //      - Find each cells distance from the boundary by running a BFS algorithm from the boundary in each connected component
    //      - Along the way, find the deepest cell in each connected component and push them to the list

    for (std::vector<int>& nBndCells : nBoundaryCellsByConnnectedComponent) {
        
        // We want to find the maximum depth cell within the connected component
        int nMaxDepthCell = nBndCells[0];

        // Run a BFS algorithm from the boundary:
        //      - All boundary cells are pushed to the queue with depth 0
        std::queue<int> nBfsQueue;
        for (int nCellIdx : nBndCells) {
            nCellDepthList[nCellIdx] = 0;
            nBfsQueue.push(nCellIdx);
        }

        while (!nBfsQueue.empty()) {
            int nCurrCell = nBfsQueue.front();
            nBfsQueue.pop();
            
            // Check if current cell is the deeper than the maximum depth cell found so far. If it is, we update nMaxDepthCell
            if (nCellDepthList[nCurrCell] > nCellDepthList[nMaxDepthCell])
                nMaxDepthCell = nCurrCell;

            // The BFS algorithm needs to be based on point-neghbours, meaning by using the Mesh-Graph
            // Note that while cells in different connected components are never face-neighbours, but they may be point-neighbours (neighbours in the Mesh-Graph)
            // Therefore, when pushing all point-neighbouring cells we need to check that they are in the same connected component
            for (int nNeiCell : nMeshGraph[nCurrCell]) {
                if (nCellDepthList[nNeiCell] != -1) continue;
                if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nCurrCell]) continue;

                nCellDepthList[nNeiCell] = nCellDepthList[nCurrCell] + 1;
                nBfsQueue.push(nNeiCell);
            }
        }

        // Finally, we can push the maximum depth cell we found
        nStartingCells.push_back(nMaxDepthCell);
    }
}


void Foam::hpathRenumber::hpathFinder::getLayers(const std::vector<int>& nStartingCells, std::vector<std::vector<int>>& nCellsByLayer) {

    // Separates all cells in the mesh to layers
    // In each connected component:
    //      1) The starting cell is set as 'Layer 0'
    //      2) All cells that are point-neighbours of the starting cell are set as layer 1
    //      3) All cells that are point-neighbours of a cell in layer 1 are set as layer 2
    //      4) Repeat step (3) with increasing layer indices until all cells have been found
    // In this way, cells in each component are grouped in layers by their MINIMUM point-distance to their starting cell

    // Note: In reality, the same BFS is run for all components at once. Because components are connected, this is equivalent 

    // Now we run BFS from starting cells
    std::queue<int> nBfsQueue;

    // Push all starting cells to queue
    for (int nCellIdx : nStartingCells) {
        nCellLayerIndex[nCellIdx] = 0;
        nBfsQueue.push(nCellIdx);
    }

    // BFS algorithm will find for each cell its minimum distance in the Mesh-Graph from the corresponding starting cell
    while(!nBfsQueue.empty())
    {
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();
        
        int nLayer = nCellLayerIndex[nCurrCell];
        if (int(nCellsByLayer.size()) <= nCellLayerIndex[nCurrCell])
            nCellsByLayer.emplace_back();
        nCellsByLayer[nLayer].push_back(nCurrCell);

        for(int nNeiCell : nMeshGraph[nCurrCell]) {
            // If neighbour is in a different connected component, ignore it
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nCurrCell]) continue;
            // If neighbour hasn't been visited yet, set it's layer and push it:
            if (nCellLayerIndex[nNeiCell] != -1) continue;
            
            nCellLayerIndex[nNeiCell] = nCellLayerIndex[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
}

void Foam::hpathRenumber::hpathFinder::solveLayer(std::vector<int> nCellsInLayer, Foam::labelList& cellOrder) {

    // Note: this method assumes all cells in the nCellsInLayer have the same 'layer index' in nCellLayerIndex
    
    // General rundown of the algorithm:
    //      1) Cells are separated into connected components
    //      2) Each connected component is solved separately and its path is added to the renumbering
    //      3) If there are still cells that have not been found, return to step (1)

    while (!nCellsInLayer.empty())
    {
        // Step 1: separate the cells into connected components
        //      - Connected components need to be connected by FACES (not points)
        std::vector<std::vector<int>> nCellsByConnnectedComponent;
        getConnectedComponents(nCellsInLayer, nCellsByConnnectedComponent);

        int nOrigFoundCellCount = nFoundCellCount;

        // Step 2: For each connected component we call getHpathinConnnectedComponent()
        for (std::vector<int>& nCellsInConnectedComponent : nCellsByConnnectedComponent)
        {
            solveConnectedComponent(nCellsInConnectedComponent, cellOrder);
        }
        
        // Find how many cells were still not found
        int nRemainingCellCount = nCellsInLayer.size() - (nFoundCellCount - nOrigFoundCellCount);
        if (nRemainingCellCount == 0) break;

        // Step 3: Find all the cells in the layer that were not found yet
        std::vector<int> nRemainingCells(nRemainingCellCount);
        int nIndex = 0;
        for (int nCellIdx : nCellsInLayer) {
            if (!bIsRenumbered[nCellIdx]) {
                nRemainingCells[nIndex++] = nCellIdx;
            }
        }

        resetCells(nRemainingCells);
        nCellsInLayer = nRemainingCells;
    }
}


void Foam::hpathRenumber::hpathFinder::getConnectedComponents(const std::vector<int>& nCellList, std::vector<std::vector<int>>& nCellsByConnnectedComponent) {
    
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
    int nConnnectedComponentIndex = -1;

    for (int nCellIdx : nCellList) {
        if (nCellConnnectedComponentIndex[nCellIdx] != -1) continue;

        // This cell hasn't been found yet
        // So we set it to be the beginning of a new connected component
        // We increment nConnnectedComponentIndex, it is now the index of this new connected component
        nConnnectedComponentIndex++;
        nCellConnnectedComponentIndex[nCellIdx] = nConnnectedComponentIndex;
        nCellsByConnnectedComponent.emplace_back();

        std::stack<int> nDfsStack;
        nDfsStack.push(nCellIdx);
        
        while(!nDfsStack.empty()) {
            int nCurrCell = nDfsStack.top();
            nDfsStack.pop();

            nCellsByConnnectedComponent[nConnnectedComponentIndex].push_back(nCurrCell);

            // We push all face-neighbouring cells that:
            //      1) Are in the same layer
            //      2) Haven't been found yet
            for (int nFaceIdx : mesh.cells()[nCurrCell]) {
                int nNeiCell = getNei(nCurrCell, nFaceIdx);
                if (nNeiCell < 0) continue;
                if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nCurrCell]) continue;
                if (nCellConnnectedComponentIndex[nNeiCell] != -1) continue;

                nCellConnnectedComponentIndex[nNeiCell] = nConnnectedComponentIndex;
                nDfsStack.push(nNeiCell);
            }
        }
    }
}


void Foam::hpathRenumber::hpathFinder::solveConnectedComponent(const std::vector<int>& nCellsInConnectedComponent, Foam::labelList& cellOrder) {

    // For small cases, find path manually (1-2 cells)
    if (nCellsInConnectedComponent.size() <= 2) {
        for (int nCellIdx : nCellsInConnectedComponent) {
            cellOrder[nFoundCellCount++] = nCellIdx;
            bIsRenumbered[nCellIdx] = true;
        }
        return;
    }
    
    // Get a graph representing the connected component
    getConnnectedComponentGraph(nCellsInConnectedComponent);

    // Get the starting cell
    int nStartCell = getStartingCellInConnnectedComponent(nCellsInConnectedComponent);

    // Reorder each cells neighbours in the ConnnectedComponent-Graph in descending order by distance from start cell
    reorderDistFromStart(nStartCell, nCellsInConnectedComponent);

    // Get a path through the connected component, using the ConnnectedComponent-Graph
    findPath(nStartCell, cellOrder);
}


void Foam::hpathRenumber::hpathFinder::getConnnectedComponentGraph(const std::vector<int>& nCellsInConnectedComponent) {

    // This method computes the ConnnectedComponent-Graph
    // Cells are connected in the ConnnectedComponent-Graph if:
    //      a) They are in the same layer
    //      b) They are face-connected in the mesh (different from Mesh-Graph - there it was point-connected)
    //          - note that this implies that they are also in the same connected component (hence the name)

    for (int nCellIdx : nCellsInConnectedComponent)
    {
        // Necessary if this being called recursively
        nConnnectedComponentGraph[nCellIdx].clear();

        // For every face on the cell, find its neighbour through that face (if there is one)
        // If that neighbour:
        //      1) Hasn't been added to the reordering yet
        //      2) Is in the same layer
        //      3) Is in the same connected component
        // Then we add it as a neighbour in the ConnectedComponent-Graph
        for (int nFaceIdx : mesh.cells()[nCellIdx])
        {
            int nNeiCell = getNei(nCellIdx, nFaceIdx);
            if (nNeiCell < 0) continue;
            if (bIsRenumbered[nNeiCell]) continue;
            if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nCellIdx]) continue;

            nConnnectedComponentGraph[nCellIdx].push_back(nNeiCell);
        }
    }
}


int Foam::hpathRenumber::hpathFinder::getStartingCellInConnnectedComponent(const std::vector<int>& nCellsInConnectedComponent) {

    // Strategy for choosing a starting cell in the connected component:
    //      1) Start from an arbitrary cell within the connected component
    //      2) Find the furthest cell from it: choose it as the starting cell

    // The reasoning behind this strategy is as follows:
    //      - In many cases the cells given may be of the general shape of a long straight path
    //      - In these cases, by starting from the farthest cell we guarantee it will be at one of the edges of the path

    // It does not matter which cell we start from, so arbitrarly start from first cell in the connected component
    int nCellIdx = nCellsInConnectedComponent[0];

    // Find farthest cell from nCellIdx
    // This is done using a standard BFS algorithm
    std::queue<int> nBfsQueue;

    bBFSFoundCell[nCellIdx] = true;
    nBfsQueue.push(nCellIdx);

    int nCurrCell = -1;

    while(!nBfsQueue.empty()) {
        nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all cells that have not yet been found by the BFS
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell]) {
            if (bBFSFoundCell[nNeiCell]) continue;
            bBFSFoundCell[nNeiCell] = true;
            nBfsQueue.push(nNeiCell);
        }
    }

    // nCurrCell should now be the farthest cell from the starting cell in the connected component
    return nCurrCell;
}


void Foam::hpathRenumber::hpathFinder::reorderDistFromStart(int nStartCell, const std::vector<int>& nCellsInConnectedComponent) {

    // findPath() tends to find much better results when each cell's neighbours are ordered in a specific way based on their face/point-distance to the starting cell

    // First, this method finds for each cell in the connected component its face-distance from the start cell
    // Secondly, this method finds for each cell in the connected component its point-distance from the start cell
    // Finally, this method reorders each cells neighbours in the ConnnectedComponent-Graph accordingly

    //  - Both of these are done using a BFS algorithm: the first using face-neigbours and the second using point-neighbours

    std::queue<int> nBfsQueue;

    // Find the face-distance of all cells in the connected component from the starting cell
    nCellFaceDistFromStart[nStartCell] = 0;
    // Push starting cell to queue
    nBfsQueue.push(nStartCell);

    while (!nBfsQueue.empty()) {
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all face-neighbours that haven't already had their face distance found
        // Face-neighbours are saved in the ConnnectedComponent-Graph
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell]) {
            if (nCellFaceDistFromStart[nNeiCell] != -1) continue;

            nCellFaceDistFromStart[nNeiCell] = nCellFaceDistFromStart[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
    
    // Find the *point*-distance of all cells in the connected component from the starting cell
    nCellPointDistFromStart[nStartCell] = 0;
    // Push starting cell to queue
    nBfsQueue.push(nStartCell);

    while (!nBfsQueue.empty()) {
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all the point-neighbour that:
        //      1) Haven't already been pushed by the BFS previously
        //      2) Are in the same layer as the Starting Cell
        //      3) Are in the same connected component as the Starting Cell
        // Point neighbours are saved in the Mesh-Graph
        for (int nNeiCell : nMeshGraph[nCurrCell]) {
            if (bIsRenumbered[nNeiCell]) continue;
            if (nCellPointDistFromStart[nNeiCell] != -1) continue;
            if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nStartCell]) continue;
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nStartCell]) continue;

            nCellPointDistFromStart[nNeiCell] = nCellPointDistFromStart[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
    
    for (int nCellIdx : nCellsInConnectedComponent) {
        // Sorting of neighbours is done in two levels:
        //      1) Neigbours are sorted *descending*-order based on their *point*-distance from the starting cell
        //      2) Neigbours with the same *point*-distance are sorted in *ascending*-order based on their *face*-distance from the starting cell
        std::sort(nConnnectedComponentGraph[nCellIdx].begin(), nConnnectedComponentGraph[nCellIdx].end(), [this](int i, int j) {
            if (nCellPointDistFromStart[i] != nCellPointDistFromStart[j])
                return nCellPointDistFromStart[i] > nCellPointDistFromStart[j];
            else
                return nCellFaceDistFromStart[i] < nCellFaceDistFromStart[j];
        });
    }
}


void Foam::hpathRenumber::hpathFinder::findPath(int nStartCell, Foam::labelList& cellOrder) {

    // Tries to find the furthest cell from the starting cell within the ConnectedComponent-Graph
    // This is done using a DFS algorithm:
    //      - Run DFS from starting cell
    //      - Return path to deepest cell found by DFS

    // Run a standard DFS algorithm from starting cell
    std::stack<int> nDfsStack;
    nDfsStack.push(nStartCell);
    nDFSParentCell[nStartCell] = nStartCell;

    // Used for finding and returning the best path
    int nBestCell = nStartCell;

    while(!nDfsStack.empty()) {
        int nCurrCell = nDfsStack.top();
        nDfsStack.pop();
        
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
        for (int nNeiIdx = nConnnectedComponentGraph[nCurrCell].size() - 1; nNeiIdx >= 0; nNeiIdx--)
        {
            int nNeiCell = nConnnectedComponentGraph[nCurrCell][nNeiIdx];
            if (nDFSCellDepth[nNeiCell] == -1) {
                // Notice we only update the nCellDepth value of a cell when we POP it from the stack
                // This means some cells may be pushed many times, but they will only be popped once
                nDFSParentCell[nNeiCell] = nCurrCell;
                nDfsStack.push(nNeiCell);
            }
        }
    }

    int nCellIdx = nBestCell;
    int nPrevCell = -1;
    while(nCellIdx != nPrevCell) {
        cellOrder[nFoundCellCount++] = nCellIdx;
        bIsRenumbered[nCellIdx] = true;

        // The first cell in the path is always its own parent
        nPrevCell = nCellIdx;
        nCellIdx = nDFSParentCell[nCellIdx];
    }
}


void Foam::hpathRenumber::hpathFinder::resetCells(std::vector<int>& nCellList) {
    // Reset the data structures for the cells that were not found in the renumbering
    for (int nCellIdx : nCellList) {
        nCellConnnectedComponentIndex[nCellIdx] = -1;
        nCellPointDistFromStart[nCellIdx] = -1;
        nCellFaceDistFromStart[nCellIdx] = -1;
        bBFSFoundCell[nCellIdx] = false;
        nDFSCellDepth[nCellIdx] = -1;
        nDFSParentCell[nCellIdx] = -1;
    }
}

