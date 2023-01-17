/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::hpathRenumber

Description
    Hpath renumber. 

SourceFiles
    hpathRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef hpathRenumber_H
#define hpathRenumber_H

#include "renumberMethod.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hpathRenumber Declaration
\*---------------------------------------------------------------------------*/

class hpathRenumber
:
    public renumberMethod
{
    // Flag indicating if layer separation is applied (default: true)
    bool m_bApplyLayerSeparation;

    //- No copy construct
    hpathRenumber(const hpathRenumber&) = delete;

    //- No copy assignment
    void operator=(const hpathRenumber&) = delete;

public:
    //- Runtime type information
    TypeName("hpath");

    // Constructors

    //- Construct given the renumber dictionary
    hpathRenumber(const dictionary& renumberDict);

    //- Destructor
    virtual ~hpathRenumber() = default;


    // Member Functions

    // Unimplemented renumber method - needs to be included for class to compile
    virtual labelList renumber(const labelListList& cellCells, const pointField& cc) const {
        NotImplemented;
        return labelList();
    };

    //  Implemented renumber method
    //  Returns the order in which cells need to be visited,
    //  from ordered back to original cell label.
    //  Uses the mesh connectivity
    virtual labelList renumber
    (
        const polyMesh& mesh,
        const pointField& cc
    ) const;


protected:

    /*---------------------------------------------------------------------------*\
                        Struct dynamicMarker Declaration
    \*---------------------------------------------------------------------------*/

    // This class supports marking and unmarking small sets of indices quickly
    // Used for efficient construction of Mesh Graph in getMeshGraph()

    class dynamicMarker {
    private:
        std::vector<bool> bMarker;
        std::stack<int> nMarkedStack;
        
    public:
        // Create a new dynamicMarker of size n
        // All items will begin unmarked
        dynamicMarker(int n);

        // Mark the i'th element (0 <= i < n), return false if it was already marked.
        bool mark(int i);

        // Clear ALL marked elements
        void clear();
    };

    /*---------------------------------------------------------------------------*\
                        Class hpathFinder Declaration
    \*---------------------------------------------------------------------------*/

    class hpathFinder {
    private: // Private member variables

        // The input mesh
        const Foam::polyMesh& mesh;
        const int nCellCount;

        // Counter for the number of renumbered cells
        int nFoundCellCount;

        // Marks cells that have been added to the renumbering as 'true'
        std::vector<bool> bIsRenumbered;

        // For every data structure, I explain what it is used for and which method is used to compute it

        // For every point a list of cells it is part of                - getPointCellLists()
        std::vector<std::vector<int>> nPntCellList;
        // For every cell a list of points on it                        - getPointCellLists()
        std::vector<std::vector<int>> nCellPntList;

        // For every cell, a list of all its point-neighbouring cells   - getMeshGraph()
        std::vector<std::vector<int>> nMeshGraph;

        // For each cell its 'layer index':                             - getLayers()  
        //      - cells in the same layer will have the same layer index
        std::vector<int> nCellLayerIndex;
        
        // For each cell its 'connected component index':                              - getConnectedComponents()
        //      - cells in the same connected component will have the same connected component index
        std::vector<int> nCellConnnectedComponentIndex;

        // For each cell its face-neighbours in the connected component                - getConnnectedComponentGraph()
        std::vector<std::vector<int>> nConnnectedComponentGraph;

        // Marks cells that have already been by BFS                                   - getStartingCellInConnnectedComponent()
        std::vector<bool> bBFSFoundCell;

        // For each cell its point distance from the connected components start cell   - reorderDistFromStart()
        std::vector<int> nCellPointDistFromStart;
        std::vector<int> nCellFaceDistFromStart;

        // For each cell its DFS depth within the connected component     - findPath()
        std::vector<int> nDFSCellDepth;
        // For each cell its DFS parent within the connected component    - findPath()
        std::vector<int> nDFSParentCell;

    public: // Public methods

        // Constructor
        hpathFinder(const Foam::polyMesh& mesh);
        
        // Get Renumbering for the mesh
        void getRenumbering(Foam::labelList& cellOrder, bool bApplyLayerSeparation);
        
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
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
