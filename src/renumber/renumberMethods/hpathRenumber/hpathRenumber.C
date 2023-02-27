/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Alon Zameret
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

#include <queue>
#include <stack>
#include <iomanip>
#include <numeric>

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

namespace
{

// This class supports marking and unmarking small sets of indices quickly
// Used for efficient construction of Mesh Graph in getMeshGraph()

class dynamicMarker
{
    std::vector<bool> bMarker;
    std::stack<int> nMarkedStack;

public:

    // Create a new dynamicMarker of size n.
    // All items will begin unmarked
    dynamicMarker(int n)
    {
        bMarker.assign(n, false);
    }

    // Mark the i'th element (0 <= i < n),
    // return false if it was already marked.
    bool mark(int i)
    {
        if (bMarker[i]) return false;
        bMarker[i] = true;
        nMarkedStack.push(i);
        return true;
    }

    // Clear ALL marked elements
    // When we want to clear all marks, we unmark every item in the stack
    // This takes amortized-time O(1), because we can only pop as many items
    // as we have pushed
    void clear()
    {
        while (!nMarkedStack.empty())
        {
            bMarker[nMarkedStack.top()] = false;
            nMarkedStack.pop();
        }
    }
};

} // End anonymous namespace


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
        const int nCellCount;

        // Counter for the number of renumbered cells
        int nFoundCellCount;

        // Marks cells that have been added to the renumbering as 'true'
        std::vector<bool> bIsRenumbered;

        // For every data structure, I explain what it is used for and which
        // method is used to compute it

        // For every point a list of cells it is part of - getPointCellLists()
        std::vector<std::vector<int>> nPntCellList;

        // For every cell a list of points on it - getPointCellLists()
        std::vector<std::vector<int>> nCellPntList;

        // For every cell, a list of all its point-neighbouring cells - getMeshGraph()
        std::vector<std::vector<int>> nMeshGraph;

        // For each cell its 'layer index': - getLayers()
        //      - cells in the same layer will have the same layer index
        std::vector<int> nCellLayerIndex;

        // For each cell its 'connected component index': - getConnectedComponents()
        //      - cells in the same connected component will have the same connected component index
        std::vector<int> nCellConnnectedComponentIndex;

        // For each cell its face-neighbours in the connected component - getConnnectedComponentGraph()
        std::vector<std::vector<int>> nConnnectedComponentGraph;

        // Marks cells that have already been by BFS - getStartingCellInConnnectedComponent()
        std::vector<bool> bBFSFoundCell;

        // For each cell its point distance from the connected components start cell   - reorderDistFromStart()
        std::vector<int> nCellPointDistFromStart;
        std::vector<int> nCellFaceDistFromStart;

        // For each cell its DFS depth within the connected component     - findPath()
        std::vector<int> nDFSCellDepth;

        // For each cell its DFS parent within the connected component    - findPath()
        std::vector<int> nDFSParentCell;


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
        int getNei(int nCellIdx, int nFaceIdx) const;

        // Creates the 'Mesh-Graph': for every cell, a list of cells that are point-neighbours with it in the mesh
        //  - Cells are point-neighbours if they have a common point
        void getMeshGraph();

        // Finds:
        //  - for each cell in the mesh a list of all points on it
        //  - for every point in the mesh a list of all cells it is part of
        void getPointCellLists();

        // Separates the mesh into layers: each cell has its layer saved in nCellLayerIndex
        // Also returns for every layer a list of all cells in it
        void getLayerSeparation(std::vector<std::vector<int>>& nCellsByLayer);

        // For every connected component, find the deepest cell and choose it as a starting cell
        // Returns a list with one starting cell per connected component
        void getStartingCells(std::vector<int>& nStartingCells) const;

        // Once the starting cells have been found, this method does the actual layer separation
        void getLayers(const std::vector<int>& nStartingCellList, std::vector<std::vector<int>>& nCellsByLayer);

        // Renumber all cells in a layer
        void solveLayer(std::vector<int> nCellsInLayer, Foam::labelList& cellOrder);

        // Seperates cells into face-connected components
        // This is done using a general DFS algorithm
        void getConnectedComponents(const std::vector<int>& nCellList, std::vector<std::vector<int>>& nCellsByConnnectedComponent);

        // Find an approximate H-path through a connected component
        void solveConnectedComponent(const std::vector<int>& nCellsInConnectedComponent, Foam::labelList& cellOrder);

        // Finds for each cell in the connected component a list of its face-neighbours within the component
        void getConnnectedComponentGraph(const std::vector<int>& nCellsInConnectedComponent);

        // Finds a starting cell within the connected component
        int getStartingCellInConnnectedComponent(const std::vector<int>& nCellsInConnectedComponent);

        // This method reorders the cells in the given connected component based on their distance from nStartCell
        void reorderDistFromStart(int nStartCell, const std::vector<int>& nCellsInConnectedComponent);

        // Finds an H-path within the connected component and returns it in nResultHpath
        // H-path is guaranteed to start at nStartCell
        void findPath(int nStartCell, Foam::labelList& cellOrder);

        // Resets data structures for the cells that weren't found
        void resetCells(std::vector<int>& nCellList);
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
    mesh(mesh),
    nCellCount(mesh.nCells())
{}


void Foam::hpathRenumber::hpathFinder::getRenumbering
(
    Foam::labelList& cellOrder,
    bool bApplyLayerSeparation
)
{
    // Find a renumbering for the entire mesh

    cellOrder.resize(nCellCount);
    // Counter for how many cells we have added to the renumbering so far
    nFoundCellCount = 0;

    // Initialize the data structures:
<<<<<<< HEAD
    std::cout << "Initializing Data Structures" << std::endl;
    initialize();
    
=======
    Info<< "Initializing Data Structures" << endl;
    initialize();

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Compute a graph to represent the mesh:
    //      - Cells in the mesh will be connected in the graph if they have a *common point*
    getMeshGraph();

    std::vector<std::vector<int>> nCellsByLayer;
<<<<<<< HEAD
    if (bApplyLayerSeparation) {
        getLayerSeparation(nCellsByLayer);
    }
    else {
=======
    if (bApplyLayerSeparation)
    {
        getLayerSeparation(nCellsByLayer);
    }
    else
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
float Foam::hpathRenumber::hpathFinder::getAccuracy(const Foam::labelList& cellOrder) const {
=======
float Foam::hpathRenumber::hpathFinder::getAccuracy
(
    const Foam::labelList& cellOrder
) const
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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
<<<<<<< HEAD
            if (nNeiIdx == nNextCellIdx) {
=======
            if (nNeiIdx == nNextCellIdx)
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
                nHitCnt++;
                break;
            }
        }
    }
    // The accuracy is the percentage of consecutive cells that were hits
    return 100.0f * float(nHitCnt) / float(int(cellOrder.size()) - 1);
}

<<<<<<< HEAD
// Private methods

void Foam::hpathRenumber::hpathFinder::initialize() {
    
=======

// Private methods
void Foam::hpathRenumber::hpathFinder::initialize()
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
int Foam::hpathRenumber::hpathFinder::getNei(int nCell, int nFaceIdx) const {
    // If the face has no neighbor, return -1
    if(nFaceIdx >= mesh.faceNeighbour().size()) 
        return -1;
    
    // Otherwise, the face connects 2 cells: the owner and the neighbor. One of these should be nCell
=======
int Foam::hpathRenumber::hpathFinder::getNei(int nCell, int nFaceIdx) const
{
    // If the face has no neighbor, return -1
    if (nFaceIdx >= mesh.faceNeighbour().size())
    {
        return -1;
    }

    // Otherwise, the face connects 2 cells: the owner and the neighbor.
    // One of these should be nCell
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    int nOwner = mesh.faceOwner()[nFaceIdx];
    int nNei = mesh.faceNeighbour()[nFaceIdx];

    // Find which of these two cells is the input cell, return the other one
    return (nOwner == nCell) ? nNei : nOwner;
}

<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getMeshGraph() {

=======

void Foam::hpathRenumber::hpathFinder::getMeshGraph()
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // First, we need to compute:
    //      1) For every point a list of cells it is part of  - nPntCellList
    //      2) For every cell a list of points on it          - nCellPntList
    getPointCellLists();

    // The purpose of the dynamic marker is to avoid checking the same cell many times in one iteration
    dynamicMarker bCellMarker(nCellCount);
<<<<<<< HEAD
    for (int nCellIdx = 0; nCellIdx < nCellCount; nCellIdx++) {
        for (int nPntIdx : nCellPntList[nCellIdx]) {
            for (int nNeiCell : nPntCellList[nPntIdx]) {
=======
    for (int nCellIdx = 0; nCellIdx < nCellCount; nCellIdx++)
    {
        for (int nPntIdx : nCellPntList[nCellIdx])
        {
            for (int nNeiCell : nPntCellList[nPntIdx])
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getPointCellLists() {

=======
void Foam::hpathRenumber::hpathFinder::getPointCellLists()
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Find:
    //  - For each cell in the mesh a list of all points on it
    nCellPntList.assign(mesh.nCells(),  std::vector<int>());
    //  - For every point in the mesh a list of all cells it is part of
    nPntCellList.assign(mesh.nPoints(), std::vector<int>());

    // The purpose of the dynamic marker is to avoid pushing the same point many times in one iteration
    dynamicMarker bPointMarker(mesh.nPoints());
    for (int nCellIdx = 0; nCellIdx < mesh.nCells(); nCellIdx++)
    {
<<<<<<< HEAD
        for(int nFaceIdx : mesh.cells()[nCellIdx]) {
            for(int nPntIdx : mesh.faces()[nFaceIdx]) {
                // If we already found nPntIdx for this cell, continue
                if (!bPointMarker.mark(nPntIdx)) continue;
                
=======
        for (int nFaceIdx : mesh.cells()[nCellIdx])
        {
            for (int nPntIdx : mesh.faces()[nFaceIdx])
            {
                // If we already found nPntIdx for this cell, continue
                if (!bPointMarker.mark(nPntIdx)) continue;

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
                // nPntIdx is on nCellIdx, so we push them to each-others lists
                nCellPntList[nCellIdx].push_back(nPntIdx);
                nPntCellList[nPntIdx].push_back(nCellIdx);
            }
        }
        // Clear the marked points for the next cell
        bPointMarker.clear();
    }
}

<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getLayerSeparation(std::vector<std::vector<int>>& nCellsByLayer)
=======
void Foam::hpathRenumber::hpathFinder::getLayerSeparation
(
    std::vector<std::vector<int>>& nCellsByLayer
)
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
{
    // We want to separate the mesh into layers:
    //      1) The mesh is separated into connected components                                          - getConnectedComponents()
    //      2) For each connected component we find the deepest cell and set it to be a starting cell   - getStartingCells()
    //      3) Each connected component is separated into layers using a BFS from the starting          - getLayers()

    std::cout << "Finding Layer Separation" << std::endl;

    // Step 1: Finding connected components
    std::vector<int> nAllCells(nCellCount);
    std::iota(nAllCells.begin(), nAllCells.end(), 0);
<<<<<<< HEAD
    
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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
<<<<<<< HEAD
    
    int nLayerCnt = nCellsByLayer.size();
    std::cout << "Mesh Separated into " << nLayerCnt << " layers" << std::endl;
    
=======

    int nLayerCnt = nCellsByLayer.size();
    std::cout << "Mesh Separated into " << nLayerCnt << " layers" << std::endl;

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Reset connected component index list for solveLayer() to use
    nCellConnnectedComponentIndex.assign(nCellCount, -1);
}

<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getStartingCells(std::vector<int>& nStartingCells) const {

    // The starting cell in each connected component should be the 'deepest' cell in that connected component
    // A cells 'depth' is its point-distance from any boundary cell
    
=======

void Foam::hpathRenumber::hpathFinder::getStartingCells
(
    std::vector<int>& nStartingCells)
const
{
    // The starting cell in each connected component should be the 'deepest' cell in that connected component
    // A cells 'depth' is its point-distance from any boundary cell

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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
<<<<<<< HEAD
        for (int nFaceIdx = range.min(); nFaceIdx <= range.max(); ++nFaceIdx) {
            for(int nPointIdx : mesh.faces()[nFaceIdx]) {
=======
        for (int nFaceIdx = range.min(); nFaceIdx <= range.max(); ++nFaceIdx)
        {
            for (int nPointIdx : mesh.faces()[nFaceIdx])
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
                bIsBoundaryPts[nPointIdx] = true;
            }
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Step 2: Finding all boundary cells
    //      - A boundary cell is a cell containing at least one boundary point
    std::vector<bool> bIsBoundaryCells(nCellCount, false);

    // For every boundary point: mark all cells that have it as boundary cells
    for (int nPntIdx = 0; nPntIdx < mesh.nPoints(); nPntIdx++)
    {
        if (!bIsBoundaryPts[nPntIdx]) continue;
        // Here we make use of the nPntCellList we found when computing the Mesh-Graph
        for (int nCellIdx : nPntCellList[nPntIdx])
<<<<<<< HEAD
            bIsBoundaryCells[nCellIdx] = true;
=======
        {
            bIsBoundaryCells[nCellIdx] = true;
        }
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    }

    // Step 3: Separate the boundary cells based on which connected component they are a part of
    //      - Find for each connected component a list of all its boundary cells
    int nConnnectedComponentCnt = *std::max_element(nCellConnnectedComponentIndex.begin(), nCellConnnectedComponentIndex.end()) + 1;
    std::vector<std::vector<int>> nBoundaryCellsByConnnectedComponent(nConnnectedComponentCnt, std::vector<int>());
    for (int nCellIdx = 0; nCellIdx < nCellCount; nCellIdx++)
    {
<<<<<<< HEAD
        if (bIsBoundaryCells[nCellIdx]) {
=======
        if (bIsBoundaryCells[nCellIdx])
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            nBoundaryCellsByConnnectedComponent[nCellConnnectedComponentIndex[nCellIdx]].push_back(nCellIdx);
        }
    }

    // Steps 4+5:
    //      - Find each cells distance from the boundary by running a BFS algorithm from the boundary in each connected component
    //      - Along the way, find the deepest cell in each connected component and push them to the list

<<<<<<< HEAD
    for (std::vector<int>& nBndCells : nBoundaryCellsByConnnectedComponent) {
        
=======
    for (std::vector<int>& nBndCells : nBoundaryCellsByConnnectedComponent)
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        // We want to find the maximum depth cell within the connected component
        int nMaxDepthCell = nBndCells[0];

        // Run a BFS algorithm from the boundary:
        //      - All boundary cells are pushed to the queue with depth 0
        std::queue<int> nBfsQueue;
<<<<<<< HEAD
        for (int nCellIdx : nBndCells) {
=======
        for (int nCellIdx : nBndCells)
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            nCellDepthList[nCellIdx] = 0;
            nBfsQueue.push(nCellIdx);
        }

<<<<<<< HEAD
        while (!nBfsQueue.empty()) {
            int nCurrCell = nBfsQueue.front();
            nBfsQueue.pop();
            
            // Check if current cell is the deeper than the maximum depth cell found so far. If it is, we update nMaxDepthCell
            if (nCellDepthList[nCurrCell] > nCellDepthList[nMaxDepthCell])
                nMaxDepthCell = nCurrCell;
=======
        while (!nBfsQueue.empty())
        {
            int nCurrCell = nBfsQueue.front();
            nBfsQueue.pop();

            // Check if current cell is the deeper than the maximum depth cell found so far. If it is, we update nMaxDepthCell
            if (nCellDepthList[nCurrCell] > nCellDepthList[nMaxDepthCell])
            {
                nMaxDepthCell = nCurrCell;
            }
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079

            // The BFS algorithm needs to be based on point-neghbours, meaning by using the Mesh-Graph
            // Note that while cells in different connected components are never face-neighbours, but they may be point-neighbours (neighbours in the Mesh-Graph)
            // Therefore, when pushing all point-neighbouring cells we need to check that they are in the same connected component
<<<<<<< HEAD
            for (int nNeiCell : nMeshGraph[nCurrCell]) {
=======
            for (int nNeiCell : nMeshGraph[nCurrCell])
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getLayers(const std::vector<int>& nStartingCells, std::vector<std::vector<int>>& nCellsByLayer) {

=======
void Foam::hpathRenumber::hpathFinder::getLayers
(
    const std::vector<int>& nStartingCells,
    std::vector<std::vector<int>>& nCellsByLayer
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Separates all cells in the mesh to layers
    // In each connected component:
    //      1) The starting cell is set as 'Layer 0'
    //      2) All cells that are point-neighbours of the starting cell are set as layer 1
    //      3) All cells that are point-neighbours of a cell in layer 1 are set as layer 2
    //      4) Repeat step (3) with increasing layer indices until all cells have been found
    // In this way, cells in each component are grouped in layers by their MINIMUM point-distance to their starting cell

<<<<<<< HEAD
    // Note: In reality, the same BFS is run for all components at once. Because components are connected, this is equivalent 
=======
    // Note: In reality, the same BFS is run for all components at once. Because components are connected, this is equivalent
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079

    // Now we run BFS from starting cells
    std::queue<int> nBfsQueue;

    // Push all starting cells to queue
<<<<<<< HEAD
    for (int nCellIdx : nStartingCells) {
=======
    for (int nCellIdx : nStartingCells)
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        nCellLayerIndex[nCellIdx] = 0;
        nBfsQueue.push(nCellIdx);
    }

    // BFS algorithm will find for each cell its minimum distance in the Mesh-Graph from the corresponding starting cell
<<<<<<< HEAD
    while(!nBfsQueue.empty())
    {
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();
        
        int nLayer = nCellLayerIndex[nCurrCell];
        if (int(nCellsByLayer.size()) <= nCellLayerIndex[nCurrCell])
            nCellsByLayer.emplace_back();
        nCellsByLayer[nLayer].push_back(nCurrCell);

        for(int nNeiCell : nMeshGraph[nCurrCell]) {
=======
    while (!nBfsQueue.empty())
    {
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        int nLayer = nCellLayerIndex[nCurrCell];
        if (int(nCellsByLayer.size()) <= nCellLayerIndex[nCurrCell])
        {
            nCellsByLayer.emplace_back();
        }
        nCellsByLayer[nLayer].push_back(nCurrCell);

        for (int nNeiCell : nMeshGraph[nCurrCell])
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            // If neighbour is in a different connected component, ignore it
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nCurrCell]) continue;
            // If neighbour hasn't been visited yet, set it's layer and push it:
            if (nCellLayerIndex[nNeiCell] != -1) continue;
<<<<<<< HEAD
            
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            nCellLayerIndex[nNeiCell] = nCellLayerIndex[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
}

<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::solveLayer(std::vector<int> nCellsInLayer, Foam::labelList& cellOrder) {

    // Note: this method assumes all cells in the nCellsInLayer have the same 'layer index' in nCellLayerIndex
    
=======

void Foam::hpathRenumber::hpathFinder::solveLayer
(
    std::vector<int> nCellsInLayer,
    Foam::labelList& cellOrder
)
{
    // Note: this method assumes all cells in the nCellsInLayer have the same 'layer index' in nCellLayerIndex

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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
<<<<<<< HEAD
        
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        // Find how many cells were still not found
        int nRemainingCellCount = nCellsInLayer.size() - (nFoundCellCount - nOrigFoundCellCount);
        if (nRemainingCellCount == 0) break;

        // Step 3: Find all the cells in the layer that were not found yet
        std::vector<int> nRemainingCells(nRemainingCellCount);
        int nIndex = 0;
<<<<<<< HEAD
        for (int nCellIdx : nCellsInLayer) {
            if (!bIsRenumbered[nCellIdx]) {
=======
        for (int nCellIdx : nCellsInLayer)
        {
            if (!bIsRenumbered[nCellIdx])
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
                nRemainingCells[nIndex++] = nCellIdx;
            }
        }

        resetCells(nRemainingCells);
        nCellsInLayer = nRemainingCells;
    }
}


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getConnectedComponents(const std::vector<int>& nCellList, std::vector<std::vector<int>>& nCellsByConnnectedComponent) {
    
=======
void Foam::hpathRenumber::hpathFinder::getConnectedComponents
(
    const std::vector<int>& nCellList,
    std::vector<std::vector<int>>& nCellsByConnnectedComponent
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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

<<<<<<< HEAD
    for (int nCellIdx : nCellList) {
=======
    for (int nCellIdx : nCellList)
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        if (nCellConnnectedComponentIndex[nCellIdx] != -1) continue;

        // This cell hasn't been found yet
        // So we set it to be the beginning of a new connected component
        // We increment nConnnectedComponentIndex, it is now the index of this new connected component
        nConnnectedComponentIndex++;
        nCellConnnectedComponentIndex[nCellIdx] = nConnnectedComponentIndex;
        nCellsByConnnectedComponent.emplace_back();

        std::stack<int> nDfsStack;
        nDfsStack.push(nCellIdx);
<<<<<<< HEAD
        
        while(!nDfsStack.empty()) {
=======

        while(!nDfsStack.empty())
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            int nCurrCell = nDfsStack.top();
            nDfsStack.pop();

            nCellsByConnnectedComponent[nConnnectedComponentIndex].push_back(nCurrCell);

            // We push all face-neighbouring cells that:
            //      1) Are in the same layer
            //      2) Haven't been found yet
<<<<<<< HEAD
            for (int nFaceIdx : mesh.cells()[nCurrCell]) {
=======
            for (int nFaceIdx : mesh.cells()[nCurrCell])
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::solveConnectedComponent(const std::vector<int>& nCellsInConnectedComponent, Foam::labelList& cellOrder) {

=======
void Foam::hpathRenumber::hpathFinder::solveConnectedComponent
(
    const std::vector<int>& nCellsInConnectedComponent,
    Foam::labelList& cellOrder
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // For small cases, find path manually (1-2 cells)
    if (nCellsInConnectedComponent.size() <= 2) {
        for (int nCellIdx : nCellsInConnectedComponent) {
            cellOrder[nFoundCellCount++] = nCellIdx;
            bIsRenumbered[nCellIdx] = true;
        }
        return;
    }
<<<<<<< HEAD
    
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Get a graph representing the connected component
    getConnnectedComponentGraph(nCellsInConnectedComponent);

    // Get the starting cell
    int nStartCell = getStartingCellInConnnectedComponent(nCellsInConnectedComponent);

    // Reorder each cells neighbours in the ConnnectedComponent-Graph in descending order by distance from start cell
    reorderDistFromStart(nStartCell, nCellsInConnectedComponent);

    // Get a path through the connected component, using the ConnnectedComponent-Graph
    findPath(nStartCell, cellOrder);
}


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::getConnnectedComponentGraph(const std::vector<int>& nCellsInConnectedComponent) {

=======
void Foam::hpathRenumber::hpathFinder::getConnnectedComponentGraph
(
    const std::vector<int>& nCellsInConnectedComponent
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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


<<<<<<< HEAD
int Foam::hpathRenumber::hpathFinder::getStartingCellInConnnectedComponent(const std::vector<int>& nCellsInConnectedComponent) {
=======
int Foam::hpathRenumber::hpathFinder::getStartingCellInConnnectedComponent
(
    const std::vector<int>& nCellsInConnectedComponent
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079

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

<<<<<<< HEAD
    while(!nBfsQueue.empty()) {
=======
    while(!nBfsQueue.empty())
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all cells that have not yet been found by the BFS
<<<<<<< HEAD
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell]) {
=======
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell])
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            if (bBFSFoundCell[nNeiCell]) continue;
            bBFSFoundCell[nNeiCell] = true;
            nBfsQueue.push(nNeiCell);
        }
    }

    // nCurrCell should now be the farthest cell from the starting cell in the connected component
    return nCurrCell;
}


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::reorderDistFromStart(int nStartCell, const std::vector<int>& nCellsInConnectedComponent) {

=======
void Foam::hpathRenumber::hpathFinder::reorderDistFromStart
(
    int nStartCell,
    const std::vector<int>& nCellsInConnectedComponent
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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

<<<<<<< HEAD
    while (!nBfsQueue.empty()) {
=======
    while (!nBfsQueue.empty())
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all face-neighbours that haven't already had their face distance found
        // Face-neighbours are saved in the ConnnectedComponent-Graph
<<<<<<< HEAD
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell]) {
=======
        for (int nNeiCell : nConnnectedComponentGraph[nCurrCell])
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            if (nCellFaceDistFromStart[nNeiCell] != -1) continue;

            nCellFaceDistFromStart[nNeiCell] = nCellFaceDistFromStart[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
<<<<<<< HEAD
    
=======

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    // Find the *point*-distance of all cells in the connected component from the starting cell
    nCellPointDistFromStart[nStartCell] = 0;
    // Push starting cell to queue
    nBfsQueue.push(nStartCell);

<<<<<<< HEAD
    while (!nBfsQueue.empty()) {
=======
    while (!nBfsQueue.empty())
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        int nCurrCell = nBfsQueue.front();
        nBfsQueue.pop();

        // Push all the point-neighbour that:
        //      1) Haven't already been pushed by the BFS previously
        //      2) Are in the same layer as the Starting Cell
        //      3) Are in the same connected component as the Starting Cell
        // Point neighbours are saved in the Mesh-Graph
<<<<<<< HEAD
        for (int nNeiCell : nMeshGraph[nCurrCell]) {
=======
        for (int nNeiCell : nMeshGraph[nCurrCell])
        {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            if (bIsRenumbered[nNeiCell]) continue;
            if (nCellPointDistFromStart[nNeiCell] != -1) continue;
            if (nCellLayerIndex[nNeiCell] != nCellLayerIndex[nStartCell]) continue;
            if (nCellConnnectedComponentIndex[nNeiCell] != nCellConnnectedComponentIndex[nStartCell]) continue;

            nCellPointDistFromStart[nNeiCell] = nCellPointDistFromStart[nCurrCell] + 1;
            nBfsQueue.push(nNeiCell);
        }
    }
<<<<<<< HEAD
    
    for (int nCellIdx : nCellsInConnectedComponent) {
        // Sorting of neighbours is done in two levels:
        //      1) Neigbours are sorted *descending*-order based on their *point*-distance from the starting cell
        //      2) Neigbours with the same *point*-distance are sorted in *ascending*-order based on their *face*-distance from the starting cell
        std::sort(nConnnectedComponentGraph[nCellIdx].begin(), nConnnectedComponentGraph[nCellIdx].end(), [this](int i, int j) {
=======

    for (int nCellIdx : nCellsInConnectedComponent)
    {
        // Sorting of neighbours is done in two levels:
        //      1) Neigbours are sorted *descending*-order based on their *point*-distance from the starting cell
        //      2) Neigbours with the same *point*-distance are sorted in *ascending*-order based on their *face*-distance from the starting cell
        std::sort
        (
            nConnnectedComponentGraph[nCellIdx].begin(),
            nConnnectedComponentGraph[nCellIdx].end(),
            [this](int i, int j) {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
            if (nCellPointDistFromStart[i] != nCellPointDistFromStart[j])
                return nCellPointDistFromStart[i] > nCellPointDistFromStart[j];
            else
                return nCellFaceDistFromStart[i] < nCellFaceDistFromStart[j];
<<<<<<< HEAD
        });
=======
            }
        );
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
    }
}


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::findPath(int nStartCell, Foam::labelList& cellOrder) {

=======
void Foam::hpathRenumber::hpathFinder::findPath
(
    int nStartCell, Foam::labelList& cellOrder
)
{
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
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

<<<<<<< HEAD
    while(!nDfsStack.empty()) {
        int nCurrCell = nDfsStack.top();
        nDfsStack.pop();
        
=======
    while (!nDfsStack.empty())
    {
        int nCurrCell = nDfsStack.top();
        nDfsStack.pop();

>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        // If the cell has already been popped previously, we can skip it
        if (nDFSCellDepth[nCurrCell] != -1) continue;

        // Otherwise, we update its depth value. This means the cell will never be PUSHED again
        nDFSCellDepth[nCurrCell] = nDFSCellDepth[nDFSParentCell[nCurrCell]] + 1;

        // If this is the new deepest cell, update best
        if (nDFSCellDepth[nCurrCell] > nDFSCellDepth[nBestCell]) {
            nBestCell = nCurrCell;
        }

<<<<<<< HEAD
        // Cells will be popped from the stack in reverse order from how we pushed them 
=======
        // Cells will be popped from the stack in reverse order from how we pushed them
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        // By iterating over the neighbors in reverse order, cells will be popped in the correct order
        for (int nNeiIdx = nConnnectedComponentGraph[nCurrCell].size() - 1; nNeiIdx >= 0; nNeiIdx--)
        {
            int nNeiCell = nConnnectedComponentGraph[nCurrCell][nNeiIdx];
<<<<<<< HEAD
            if (nDFSCellDepth[nNeiCell] == -1) {
=======
            if (nDFSCellDepth[nNeiCell] == -1)
            {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
                // Notice we only update the nCellDepth value of a cell when we POP it from the stack
                // This means some cells may be pushed many times, but they will only be popped once
                nDFSParentCell[nNeiCell] = nCurrCell;
                nDfsStack.push(nNeiCell);
            }
        }
    }

    int nCellIdx = nBestCell;
    int nPrevCell = -1;
<<<<<<< HEAD
    while(nCellIdx != nPrevCell) {
=======
    while (nCellIdx != nPrevCell)
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        cellOrder[nFoundCellCount++] = nCellIdx;
        bIsRenumbered[nCellIdx] = true;

        // The first cell in the path is always its own parent
        nPrevCell = nCellIdx;
        nCellIdx = nDFSParentCell[nCellIdx];
    }
}


<<<<<<< HEAD
void Foam::hpathRenumber::hpathFinder::resetCells(std::vector<int>& nCellList) {
    // Reset the data structures for the cells that were not found in the renumbering
    for (int nCellIdx : nCellList) {
=======
void Foam::hpathRenumber::hpathFinder::resetCells
(
    std::vector<int>& nCellList
)
{
    // Reset the data structures for the cells that were not found in the renumbering
    for (int nCellIdx : nCellList)
    {
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
        nCellConnnectedComponentIndex[nCellIdx] = -1;
        nCellPointDistFromStart[nCellIdx] = -1;
        nCellFaceDistFromStart[nCellIdx] = -1;
        bBFSFoundCell[nCellIdx] = false;
        nDFSCellDepth[nCellIdx] = -1;
        nDFSParentCell[nCellIdx] = -1;
    }
}

<<<<<<< HEAD
=======

// ************************************************************************* //
>>>>>>> 69df0ad468d9546727669151caaddbff149e8079
