/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::DistributedDelaunayMesh

Description

SourceFiles
    DistributedDelaunayMeshI.H
    DistributedDelaunayMesh.C
    DistributedDelaunayMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef DistributedDelaunayMesh_H
#define DistributedDelaunayMesh_H

#include "DelaunayMesh.H"
#include "backgroundMeshDecomposition.H"
#include "autoPtr.H"
#include "boundBox.H"
#include "indexedVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;

/*---------------------------------------------------------------------------*\
                  Class DistributedDelaunayMesh Declaration
\*---------------------------------------------------------------------------*/

template<class Triangulation>
class DistributedDelaunayMesh
:
    public DelaunayMesh<Triangulation>
{
public:

        typedef typename Triangulation::Vertex_handle Vertex_handle;
        typedef typename Triangulation::Cell_handle   Cell_handle;
        typedef typename Triangulation::Point         Point;

        typedef typename Triangulation::Finite_vertices_iterator
            Finite_vertices_iterator;
        typedef typename Triangulation::Finite_cells_iterator
            Finite_cells_iterator;
        typedef typename Triangulation::All_cells_iterator
            All_cells_iterator;

private:

        autoPtr<List<boundBox>> allBackgroundMeshBounds_;


    // Private Member Functions

        //-
        bool distributeBoundBoxes(const boundBox& bb);

        //-
        bool isLocal(const Vertex_handle& v) const;

        bool isLocal(const label localProcIndex) const;

        labelList overlapProcessors
        (
            const point& centre,
            const scalar radiusSqr
        ) const;

        bool checkProcBoundaryCell
        (
            const Cell_handle& cit,
            Map<labelList>& circumsphereOverlaps
        ) const;

        void findProcessorBoundaryCells
        (
            Map<labelList>& circumsphereOverlaps
        ) const;

        void markVerticesToRefer
        (
            const Map<labelList>& circumsphereOverlaps,
            PtrList<labelPairHashSet>& referralVertices,
            DynamicList<label>& targetProcessor,
            DynamicList<Vb>& parallelInfluenceVertices
        );

        label referVertices
        (
            const labelUList& targetProcessor,
            DynamicList<Vb>& parallelVertices,
            PtrList<labelPairHashSet>& referralVertices,
            labelPairHashSet& receivedVertices
        );

        //- No copy construct
        DistributedDelaunayMesh
        (
            const DistributedDelaunayMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const DistributedDelaunayMesh<Triangulation>&) = delete;


public:

    // Constructors

        //- Construct from components
        explicit DistributedDelaunayMesh(const Time& runTime);

        DistributedDelaunayMesh
        (
            const Time& runTime,
            const word& meshName
        );


    //- Destructor
    ~DistributedDelaunayMesh() = default;


    // Queries

        //- Use DelaunayMesh timeCheck function
        using DelaunayMesh<Triangulation>::timeCheck;

        scalar calculateLoadUnbalance() const;


    // Member Functions

        //- Build a mapDistribute for the supplied destination processor data
        static autoPtr<mapDistribute> buildMap(const labelUList& toProc);

        //-
        bool distribute(const boundBox& bb);

        autoPtr<mapDistribute> distribute
        (
            const backgroundMeshDecomposition& decomposition,
            List<Foam::point>& points
        );

        //- Refer vertices so that the processor interfaces are consistent
        void sync(const boundBox& bb);

        //- Refer vertices so that the processor interfaces are consistent
        void sync
        (
            const boundBox& bb,
            PtrList<labelPairHashSet>& referralVertices,
            labelPairHashSet& receivedVertices,
            bool iterateReferral = true
        );

        //- Inserts points into the triangulation if the point is within
        //  the circumsphere of another cell. Returns HashSet of failed
        //  point insertions
        template<class PointIterator>
        labelPairHashSet rangeInsertReferredWithInfo
        (
            PointIterator begin,
            PointIterator end,
            bool printErrors = true
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "DistributedDelaunayMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DistributedDelaunayMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
