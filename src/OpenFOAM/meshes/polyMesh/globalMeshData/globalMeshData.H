/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::globalMeshData

Description
    Various mesh related information for a parallel run. Upon construction,
    constructs all info using parallel communication.

    Requires:
    - all processor patches to have correct ordering.
    - all processorPatches to have their transforms set.

    The shared point and edge addressing calculates addressing for points
    and edges on coupled patches.  In the 'old' way a distinction was made
    between points/edges that are only on two processors and those that are
    on multiple processors.  The problem is that those on multiple processors
    do not allow any transformations and require a global reduction on the
    master processor.

    The alternative is to have an exchange schedule (through a 'mapDistribute')
    which sends all point/edge data (no distinction is made between
    those on two and those on more than two coupled patches) to the local
    'master'.  This master then does any calculation and sends
    the result back to the 'slave' points/edges.  This only needs to be done
    on points on coupled faces.  Any transformation is done using a
    predetermined set of transformations - since transformations have to be
    space filling only a certain number of transformation is supported.

    The exchange needs
    - a field of data
    - a mapDistribute which does all parallel exchange and transformations
      This appends remote data to the end of the field.
    - a set of indices which indicate where to get untransformed data in the
      field
    - a set of indices which indicate where to get transformed data in the
      field

Note
    - compared to 17x nTotalFaces, nTotalPoints do not compensate for
      shared points since this would trigger full connectivity analysis
    - most calculation is demand driven and uses parallel communication
      so make sure to invoke on all processors at the same time
    - old sharedEdge calculation: currently an edge is considered shared
      if it uses two shared points and is used more than once.  This is not
      correct on processor patches but it only slightly overestimates the number
      of shared edges.  Doing full analysis of how many patches use the edge
      would be too complicated

See also
    mapDistribute
    globalIndexAndTransform

SourceFiles
    globalMeshData.C
    globalMeshDataTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_globalMeshData_H
#define Foam_globalMeshData_H

#include "processorTopology.H"
#include "globalIndex.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class mapDistribute;
template<class T> class CompactListList;
template<class T> class EdgeMap;
class globalIndexAndTransform;
class bitSet;

/*---------------------------------------------------------------------------*\
                      Class globalMeshData Declaration
\*---------------------------------------------------------------------------*/

class globalMeshData
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Data related to the complete mesh

        //- Global numbering for mesh points.
        //- Not compensated for duplicate points!
        globalIndex globalMeshPointAddr_;

        //- Global numbering for mesh faces.
        //- Not compensated for duplicate faces!
        globalIndex globalMeshFaceAddr_;

        //- Global numbering for mesh cells.
        globalIndex globalMeshCellAddr_;


        // Processor patch addressing (be careful if not running in parallel!)

            //- The processor/processor topology
            processorTopology processorTopology_;

            //- List of processor patch labels
            //  (size of list = number of processor patches)
            labelList processorPatches_;

            //- List of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor patches.
            //  (size of list = number of patches)
            labelList processorPatchIndices_;

            //- processorPatchIndices_ of the neighbours processor patches
            labelList processorPatchNeighbours_;


        // Coupled point addressing
        // This is addressing from coupled point to coupled points/faces/cells.
        // This is a full schedule so includes points used by only two
        // coupled patches.

            //- Patch of coupled faces. Additional patch edge to mesh edges
            //  correspondence:
            //      points: meshPoints(), meshPointMap()
            //      edges : meshEdges(), meshEdgeMap()
            mutable autoPtr<indirectPrimitivePatch> coupledPatchPtr_;
            mutable autoPtr<labelList> coupledPatchMeshEdgesPtr_;
            mutable autoPtr<Map<label>> coupledPatchMeshEdgeMapPtr_;

            //- Global numbering for coupledPatch points
            mutable autoPtr<globalIndex> globalPointNumberingPtr_;

            //- Global numbering for transforms
            mutable autoPtr<globalIndexAndTransform> globalTransformsPtr_;

            // Coupled point to coupled points

            mutable autoPtr<labelListList> globalPointSlavesPtr_;
            mutable autoPtr<labelListList> globalPointTransformedSlavesPtr_;
            mutable autoPtr<mapDistribute> globalPointSlavesMapPtr_;

            // Coupled edge to coupled edges

            mutable autoPtr<globalIndex> globalEdgeNumberingPtr_;
            mutable autoPtr<labelListList> globalEdgeSlavesPtr_;
            mutable autoPtr<labelListList> globalEdgeTransformedSlavesPtr_;
            mutable autoPtr<bitSet> globalEdgeOrientationPtr_;
            mutable autoPtr<mapDistribute> globalEdgeSlavesMapPtr_;


            // Coupled point to boundary faces

            mutable autoPtr<globalIndex> globalBoundaryFaceNumberingPtr_;
            mutable autoPtr<labelListList> globalPointBoundaryFacesPtr_;
            mutable autoPtr<labelListList>
                globalPointTransformedBoundaryFacesPtr_;
            mutable autoPtr<mapDistribute> globalPointBoundaryFacesMapPtr_;

            // Coupled point to boundary cells

            mutable autoPtr<labelList> boundaryCellsPtr_;
            mutable autoPtr<globalIndex> globalBoundaryCellNumberingPtr_;
            mutable autoPtr<labelListList> globalPointBoundaryCellsPtr_;
            mutable autoPtr<labelListList>
                globalPointTransformedBoundaryCellsPtr_;
            mutable autoPtr<mapDistribute> globalPointBoundaryCellsMapPtr_;


            // Other: coupled point to coupled COLLOCATED points
            mutable autoPtr<labelListList> globalCoPointSlavesPtr_;
            mutable autoPtr<mapDistribute> globalCoPointSlavesMapPtr_;



        // Globally shared point addressing

            //- Total number of global points
            mutable label nGlobalPoints_;

            //- Indices of local points that are globally shared
            mutable autoPtr<labelList> sharedPointLabelsPtr_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            mutable autoPtr<labelList> sharedPointAddrPtr_;

            //- Shared point global labels.
            //  Global point index for every local shared point.
            //  Only valid if constructed with this information or if
            //  pointProcAddressing read.
            mutable autoPtr<labelList> sharedPointGlobalLabelsPtr_;


        // Globally shared edge addressing. Derived from shared points.
        // All demand driven since don't want to construct edges always.

            //- Total number of global edges
            mutable label nGlobalEdges_;

            //- Indices of local edges that are globally shared
            mutable autoPtr<labelList> sharedEdgeLabelsPtr_;

            //- Indices of globally shared edge in the master list
            //  This list contains all the shared edges in the mesh
            mutable autoPtr<labelList> sharedEdgeAddrPtr_;


    // Private Member Functions

        //- Set up processor patch addressing
        void initProcAddr();

        //- Helper function for shared edge addressing
        static void countSharedEdges
        (
            const EdgeMap<labelList>&,
            EdgeMap<label>&,
            label&
        );

        //- Calculate shared point addressing
        void calcSharedPoints() const;

        //- Calculate shared edge addressing
        void calcSharedEdges() const;

        //- Calculate global point addressing.
        void calcGlobalPointSlaves() const;

        // Global edge addressing

            //- Calculate connected points
            void calcPointConnectivity(List<labelPairList>&) const;

            //- Calculate pointEdges and pointPoints addressing
            void calcGlobalPointEdges
            (
                labelListList& globalPointEdges,
                List<labelPairList>& globalPointPoints
            ) const;

            //- Look up remote and local point and find using info the
            //  transforms to go from remotePoint to localPoint
            label findTransform
            (
                const labelPairList& info,
                const labelPair& remotePoint,
                const label localPoint
            ) const;

            //- Calculate global edge addressing.
            void calcGlobalEdgeSlaves() const;

            //- Calculate orientation w.r.t. edge master.
            void calcGlobalEdgeOrientation() const;


        // Global boundary face/cell addressing

            //- Calculate coupled point to uncoupled boundary faces. Local only.
            void calcPointBoundaryFaces(labelListList&) const;

            //- Calculate global point to global boundary face addressing.
            void calcGlobalPointBoundaryFaces() const;

            //- Calculate global point to global boundary cell addressing.
            void calcGlobalPointBoundaryCells() const;

        // Other

            // Point to collocated points. Note that not all points on
            // coupled patches now have a master! (since points on either
            // side of a cyclic are not connected). So check whether the map
            // reaches all points and decide who is master, slave and who is
            // its own master. Maybe store as well?

            void calcGlobalCoPointSlaves() const;

public:

    // Generated Methods

        //- No copy construct
        globalMeshData(const globalMeshData&) = delete;

        //- No copy assignment
        void operator=(const globalMeshData&) = delete;


    //- Runtime type information
    ClassName("globalMeshData");


    // Static data members

        //- Geometric tolerance (fraction of bounding box)
        static const Foam::scalar matchTol_;


    // Constructors

        //- Construct from mesh, derive rest (does parallel communication!)
        globalMeshData(const polyMesh& mesh);


    //- Destructor
    ~globalMeshData();

        //- Remove all demand driven data
        void clearOut();


    // Member Functions

        // Access

            //- Return the mesh reference
            const polyMesh& mesh() const noexcept
            {
                return mesh_;
            }

            //- Does the mesh contain processor patches?
            //- (also valid when not running parallel)
            bool parallel() const noexcept
            {
                return !processorPatches_.empty();
            }

            //- Global numbering for \b mesh points.
            //- Not compensated for duplicate points!
            const globalIndex& globalMeshPointAddr() const noexcept
            {
                return globalMeshPointAddr_;
            }

            //- Global numbering for \b mesh faces.
            //- Not compensated for duplicate faces!
            const globalIndex& globalMeshFaceAddr() const noexcept
            {
                return globalMeshFaceAddr_;
            }

            //- Global numbering for \b mesh cells.
            const globalIndex& globalMeshCellAddr() const noexcept
            {
                return globalMeshCellAddr_;
            }

            //- Total global number of \b mesh points.
            //- Not compensated for duplicate points!
            label nTotalPoints() const noexcept
            {
                return globalMeshPointAddr_.totalSize();
            }

            //- Total global number of \b mesh faces.
            //- Not compensated for duplicate faces!
            label nTotalFaces() const noexcept
            {
                return globalMeshFaceAddr_.totalSize();
            }

            //- Total global number of \b mesh cells.
            label nTotalCells() const noexcept
            {
                return globalMeshCellAddr_.totalSize();
            }


        // Processor patch addressing (be careful when not running in parallel)

            //- The processor to processor topology.
            const processorTopology& topology() const noexcept
            {
                return processorTopology_;
            }

            //- Order in which the patches should be initialised/evaluated
            //- corresponding to the schedule
            const lduSchedule& patchSchedule() const noexcept
            {
                return processorTopology_.patchSchedule();
            }

            //- Return list of processor patch labels
            //  (size of list = number of processor patches)
            const labelList& processorPatches() const noexcept
            {
                return processorPatches_;
            }

            //- Return list of indices into processorPatches_ for each patch.
            //  Index = -1 for non-processor parches.
            //  (size of list = number of patches)
            const labelList& processorPatchIndices() const noexcept
            {
                return processorPatchIndices_;
            }

            //- Return processorPatchIndices of the neighbours
            //- processor patches. -1 if not running parallel.
            const labelList& processorPatchNeighbours() const noexcept
            {
                return processorPatchNeighbours_;
            }


        // Globally shared point addressing

            //- Return number of globally shared points
            label nGlobalPoints() const;

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const;

            //- Return addressing into the complete globally shared points
            //- list
            //  Note: It is assumed that a (never constructed) complete
            //  list of globally shared points exists.  The set of shared
            //  points on the current processor is a subset of all shared
            //  points. Shared point addressing gives the index in the
            //  list of all globally shared points for each of the locally
            //  shared points.
            const labelList& sharedPointAddr() const;

            //- Return shared point global labels. Tries to read
            //- 'pointProcAddressing' and returns list or -1 if none
            //- available.
            const labelList& sharedPointGlobalLabels() const;

            //- Collect coordinates of shared points on all processors.
            //- (does parallel communication!)
            //  Note: not valid for cyclicParallel since shared cyclic points
            //  are merged into single global point. (use geometricSharedPoints
            //  instead)
            pointField sharedPoints() const;

            //- Like sharedPoints but keeps cyclic points separate.
            //- (does geometric merging; uses matchTol_*bb as merging tolerance)
            //- Use sharedPoints() instead.
            pointField geometricSharedPoints() const;



        // Globally shared edge addressing

            //- Return number of globally shared edges.
            //  Demand-driven
            //  calculation so call needs to be synchronous among processors!
            label nGlobalEdges() const;

            //- Return indices of local edges that are globally shared.
            //  Demand-driven
            //  calculation so call needs to be synchronous among processors!
            const labelList& sharedEdgeLabels() const;

            //- Return addressing into the complete globally shared edge list.
            //  The set of shared
            //  edges on the current processor is a subset of all shared
            //  edges. Shared edge addressing gives the index in the
            //  list of all globally shared edges for each of the locally
            //  shared edges.
            //  Demand-driven
            //  calculation so call needs to be synchronous among processors!
            const labelList& sharedEdgeAddr() const;



        // Global master - slave point communication

            //- Return patch of all coupled faces
            const indirectPrimitivePatch& coupledPatch() const;

            //- Return map from coupledPatch edges to mesh edges
            const labelList& coupledPatchMeshEdges() const;

            //- Return map from mesh edges to coupledPatch edges
            const Map<label>& coupledPatchMeshEdgeMap() const;

            //- Global transforms numbering
            const globalIndexAndTransform& globalTransforms() const;

            //- Helper: synchronise data with transforms
            template<class Type, class CombineOp, class TransformOp>
            static void syncData
            (
                List<Type>& elems,
                const labelListList& slaves,
                const labelListList& transformedSlaves,
                const mapDistribute& slavesMap,
                const globalIndexAndTransform&,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Helper: synchronise data without transforms
            template<class Type, class CombineOp>
            static void syncData
            (
                List<Type>& elems,
                const labelListList& slaves,
                const labelListList& transformedSlaves,
                const mapDistribute& slavesMap,
                const CombineOp& cop
            );


            // Coupled point to coupled points. Coupled points are
            // points on any coupled patch.

                //- Numbering of coupled points is according to coupledPatch.
                const globalIndex& globalPointNumbering() const;
                const labelListList& globalPointSlaves() const;
                const labelListList& globalPointTransformedSlaves() const;
                const mapDistribute& globalPointSlavesMap() const;
                //- Helper to synchronise coupled patch point data
                template<class Type, class CombineOp, class TransformOp>
                void syncPointData
                (
                    List<Type>& pointData,
                    const CombineOp& cop,
                    const TransformOp& top
                ) const;

            // Coupled edge to coupled edges.

                const globalIndex& globalEdgeNumbering() const;
                const labelListList& globalEdgeSlaves() const;
                const labelListList& globalEdgeTransformedSlaves() const;
                const mapDistribute& globalEdgeSlavesMap() const;
                //- Is my edge same orientation as master edge
                const bitSet& globalEdgeOrientation() const;

            // Collocated point to collocated point

                const labelListList& globalCoPointSlaves() const;
                const mapDistribute& globalCoPointSlavesMap() const;

            // Coupled point to boundary faces. These are uncoupled boundary
            // faces only but include empty patches.

                //- Numbering of boundary faces is face-mesh.nInternalFaces()
                const globalIndex& globalBoundaryFaceNumbering() const;
                const labelListList& globalPointBoundaryFaces() const;
                const labelListList& globalPointTransformedBoundaryFaces()
                const;
                const mapDistribute& globalPointBoundaryFacesMap() const;

            // Coupled point to boundary cell

                //- From boundary cell to mesh cell
                const labelList& boundaryCells() const;

                //- Numbering of boundary cells is according to boundaryCells()
                const globalIndex& globalBoundaryCellNumbering() const;
                const labelListList& globalPointBoundaryCells() const;
                const labelListList& globalPointTransformedBoundaryCells()
                const;
                const mapDistribute& globalPointBoundaryCellsMap() const;


            // Other

                //- Helper for merging (collocated!) mesh point data.
                //  Determines:
                //  - my unique indices
                //  - global numbering over all unique indices
                //  - the global number for all local points (so this will
                //    be local for my unique points)
                autoPtr<globalIndex> mergePoints
                (
                    labelList& pointToGlobal,
                    labelList& uniquePoints
                ) const;

                //- Helper for merging (collocated!) patch point data.
                //  Takes maps from:
                //  local points to/from mesh. Determines
                //  - my unique points. These are mesh point indices, not patch
                //    point indices.
                //  - global numbering over all unique indices.
                //  - the global number for all local points.
                autoPtr<globalIndex> mergePoints
                (
                    const labelUList& meshPoints,
                    const Map<label>& meshPointMap,  //!< currently unused
                    labelList& pointToGlobal,
                    labelList& uniqueMeshPoints
                ) const;


        // Edit

            //- Update for moving points.
            void movePoints(const pointField& newPoints);

            //- Change global mesh data given a topological change. Does a
            //  full parallel analysis to determine shared points and
            //  boundaries.
            void updateMesh();


    // Mesh Topology Calculation

        //- Determine (local or global) cellCells from mesh agglomeration.
        //  Agglomeration is local to the processor.
        //
        //  - parallel = false
        //    Resulting connections are in local cell indices.
        //    Coupled across cyclics but not processor patches.
        //  - parallel = true
        //    Resulting connections are in global cell indices.
        //    Coupled across cyclics and processor patches.
        static void calcCellCells
        (
            const polyMesh& mesh,
            //! The cell agglomeration, negative agglomeration are excluded
            //! from the sub-mesh
            const labelUList& agglom,
            //! The number of unique, agglomerated cells
            const label nLocalCoarse,
            //! Use global cell ids in parallel
            const bool parallel,
            //! [out] the mesh connectivity as CSR
            CompactListList<label>& cellCells
        );

        //- Determine (local or global) cellCells and face weights
        //- from mesh agglomeration.
        //  Uses mag of faceArea as weights
        static void calcCellCells
        (
            const polyMesh& mesh,
            //! The cell agglomeration, negative agglomeration are excluded
            //! from the sub-mesh
            const labelUList& agglom,
            //! The number of unique, agglomerated cells
            const label nLocalCoarse,
            //! Use global cell ids in parallel
            const bool parallel,
            //! [out] the mesh connectivity as CSR
            CompactListList<label>& cellCells,
            //! [out] the connectivity weights (face area)
            CompactListList<scalar>& cellCellWeights
        );

        //- Determine (local or global) mesh connectivity
        //
        //  - parallel = false
        //    Resulting connections are in local cell indices.
        //    Coupled across cyclics but not processor patches.
        //  - parallel = true
        //    Resulting connections are in global cell indices.
        //    Coupled across cyclics and processor patches.
        static void calcCellCells
        (
            const polyMesh& mesh,
            //! [out] the mesh connectivity as CSR
            CompactListList<label>& cellCells,
            //! Use global cell ids in parallel
            const bool parallel = false
        );

        //- Determine (local or global) sub-mesh connectivity
        //
        //  - parallel = false
        //    Resulting connections are in local cell indices.
        //    Coupled across cyclics but not processor patches.
        //  - parallel = true
        //    Resulting connections are in global cell indices.
        //    Coupled across cyclics and processor patches.
        //
        //  \return the cellMap
        static labelList calcCellCells
        (
            const polyMesh& mesh,
            //! The cell-subset
            const bitSet& selectedCells,
            //! [out] the mesh connectivity as CSR
            CompactListList<label>& cellCells,
            //! Use global cell ids in parallel
            const bool parallel = false
        );

        //- Determine (local or global) sub-mesh connectivity
        //
        //  - parallel = false
        //    Resulting connections are in local cell indices.
        //    Coupled across cyclics but not processor patches.
        //  - parallel = true
        //    Resulting connections are in global cell indices.
        //    Coupled across cyclics and processor patches.
        //
        //  \return the cellMap
        static labelList calcCellCells
        (
            const polyMesh& mesh,
            //! The cell-subset
            const labelUList& selectedCells,
            //! [out] the mesh connectivity as CSR
            CompactListList<label>& cellCells,
            //! Use global cell ids in parallel
            const bool parallel = false
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalMeshDataTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
