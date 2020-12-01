/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::distanceSurface

Description
    A surface defined by a distance from an input searchable surface.
    Uses an iso-surface algorithm (cell, topo, point) for constructing the
    distance surface.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type        distanceSurface;
            distance    0;

            surfaceType triSurfaceMesh;
            surfaceName something.obj;

            topology    nearestPoints;
            nearestPoints
            (
                (0 0 0)
                (10 10 0)
            );

            // Max distance from nearestPoints
            maxDistance 0.005;
        }
    }
    \endverbatim

    Dictionary controls:
    \table
        Property | Description                              | Required | Default
        distance | distance from surface                    | no  | 0
        signed   | Use sign when distance is positive       | no  | true
        isoAlgorithm | (cell/point/topo)                    | no  | topo
        regularise | Face simplification (enum or bool)     | no  | true
        bounds   | Limit with bounding box                  | no  |
        surfaceType | Type of surface                       | yes |
        surfaceName | Name of surface in \c triSurface/     | no  | dict name
        topology    | Topology filter (none/largestRegion/nearestPoints/proximity) | no  | -
        nearestPoints | Point-based surface segmentation    | no  | -
        maxDistance | Max search distance for nearestPoints | no  | GREAT
        relProximity | Max limit of absolute vs normal distance | no  | 10
    \endtable

    Filtering

    - \c largestRegion (pre-filter):
        The cut cells are checked for topological connectivity and the
        region with the most number of cut cells is retained.
        This handles the "ragged" edge problem.

    - \c nearestPoints (pre-filter):
        The cut cells split into regions, the regions closest to the
        user-defined points are retained.
        There is no max search

    - \c proximity (post-filter):
        Checks the resulting faces against the original search surface
        and rejects faces with an absolute distance that is greater
        than (relProximity * normal distance).
        This traps the ragged overhanging cut faces that do not correspond
        to the originating surface normals.
    .

Note
    For distance = 0, some special adjustments.
    - Always signed (ignoring the input value).
    - Use normal distance from surface (for better treatment of open edges).
    - Additional checks for open surfaces edges are used to limit the extend
      of resulting distance surface.
      The resulting surface elements will, however, contain partial cell
      coverage. NB: Not applicable if the \c point isoAlgorithm is used.

The keyword \c cell (bool value) which was use in 1906 and earlier to switch
between point/cell algorithms is now ignored (2020-12).

In (2020-12) changed the default algorithm from cell to topo.

SourceFiles
    distanceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef distanceSurface_H
#define distanceSurface_H

#include "sampledSurface.H"
#include "searchableSurface.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class distanceSurface Declaration
\*---------------------------------------------------------------------------*/

class distanceSurface
{
    // Data Types

        //- The type of cell/face cuts
        enum class topologyFilterType : uint8_t
        {
            NONE,               //!< No additional filtering
            LARGEST_REGION,     //!< Retain largest region
            NEAREST_POINTS,     //!< Retain regions nearest to the points
            PROXIMITY           //!< Post-filter for surface proximity
        };

        //- Names for the topology filter
        static const Enum<topologyFilterType> topoFilterNames_;


    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Searchable surface
        const autoPtr<searchableSurface> geometryPtr_;

        //- Distance value
        const scalar distance_;

        //- Distance is zero. Implies signed and additional optimizations
        const bool withZeroDistance_;

        //- Use signed distance
        const bool withSignDistance_;

        //- Parameters for iso-surface (algorithm, filter, mergeTol, etc)
        isoSurfaceParams isoParams_;

        //- Optional topology filter
        topologyFilterType topoFilter_;

        //- Points for nearest-points segmentation
        pointField nearestPoints_;

        //- Max search distance squared (for nearestPoints)
        scalar maxDistanceSqr_;

        //- Max search distance squared (for nearestPoints)
        scalar relProximity_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;


    // Sampling geometry. Directly stored or via an iso-surface (ALGO_POINT)

        //- The extracted surface (direct storage)
        mutable meshedSurface surface_;

        //- For every face the original cell in mesh (direct storage)
        mutable labelList meshCells_;

        //- Extracted iso-surface, for interpolators
        mutable autoPtr<isoSurfaceBase> isoSurfacePtr_;


protected:

    // Protected Member Functions

        //- Is currently backed by an isoSurfacePtr_
        bool hasIsoSurface() const
        {
            return bool(isoSurfacePtr_);
        }

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> isoSurfaceInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& cellValues,
            const Field<Type>& pointValues
        ) const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->interpolate(cellValues, pointValues);
            }

            return nullptr;
        }


    // Private Member Functions

        //- Prepare for region split.
        bool refineBlockedCells
        (
            bitSet& ignoreCells,
            const isoSurfaceBase& isoContext
        ) const;

        //- Prepare for region split.
        bitSet filterPrepareRegionSplit
        (
            const bitSet& ignoreCells
        ) const;

        //- Keep region with the most cuts (after region split)
        void filterKeepLargestRegion(bitSet& ignoreCells) const;

        //- Keep region(s) closest to the nearest points
        void filterKeepNearestRegions(bitSet& ignoreCells) const;

        //- Adjust extracted iso-surface to remove far faces
        void filterByProximity();


public:

    //- Runtime type information
    TypeName("distanceSurface");


    // Constructors

        //- Construct from dictionary
        distanceSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from components with zero-distanced
        distanceSurface
        (
            const polyMesh& mesh,
            const word& surfaceType,
            const word& surfaceName,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bool interpolate = false
        );

        //- Construct from components
        distanceSurface
        (
            const polyMesh& mesh,
            const bool interpolate,
            const word& surfaceType,
            const word& surfaceName,
            const scalar distance,
            const bool useSignedDistance,
            const isoSurfaceParams& params = isoSurfaceParams()
        );


    //- Destructor
    virtual ~distanceSurface() = default;


    // Member Functions

        //- Create/recreate the distance surface
        void createGeometry();

        //- The name of the underlying searchableSurface
        const word& surfaceName() const
        {
            return geometryPtr_->name();
        }

        //- The distance to the underlying searchableSurface
        scalar distance() const noexcept
        {
            return distance_;
        }

        //- The underlying surface
        const meshedSurface& surface() const
        {
            if (isoSurfacePtr_)
            {
                return *isoSurfacePtr_;
            }
            return surface_;
        }

        //- The underlying surface
        meshedSurface& surface()
        {
            if (isoSurfacePtr_)
            {
                return *isoSurfacePtr_;
            }
            return surface_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return meshCells_;
        }


    // Output

        //- Print information
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
