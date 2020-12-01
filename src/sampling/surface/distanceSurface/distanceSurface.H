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

            //- Max distance from nearestPoints
            maxDistance 0.005;
        }
    }
    \endverbatim

    Dictionary controls:
    \table
        Property | Description                              | Required | Default
        distance | distance from surface                    | yes |
        signed   | Use sign when distance is positive       | no  | true
        isoAlgorithm | (cell/point/topo)                    | no  | topo
        regularise | Face simplification (enum or bool)     | no  | true
        bounds   | Limit with bounding box                  | no  |
        surfaceType | Type of surface                       | yes |
        surfaceName | Name of surface in \c triSurface/     | no  | dict name
        topology    | Topology filter (none/largestRegion/nearestPoints) | no  | -
        nearestPoints | Point-based surface segmentation    | no  | -
        maxDistance | Max search distance for nearestPoints | no  | GREAT
    \endtable

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
        };

        //- Names for the topology filter
        static const Enum<topologyFilterType> topoFilterNames_;


    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Searchable surface
        const autoPtr<searchableSurface> surfacePtr_;

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

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;

        //- Constructed iso surface
        autoPtr<isoSurfaceBase> isoSurfacePtr_;


    // Private Member Functions

        //- Prepare for region split.
        bool refineBlockedCells
        (
            bitSet& ignoreCells,
            const isoSurfaceBase& isoContext
        ) const;

        //- Prepare for region split.
        boolList filterPrepareRegionSplit
        (
            const bitSet& ignoreCells
        ) const;

        //- Keep region with the most cuts (after region split)
        void filterKeepLargestRegion(bitSet& ignoreCells) const;

        //- Keep region(s) closest to the nearest points
        void filterKeepNearestRegions(bitSet& ignoreCells) const;


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
            return surfacePtr_->name();
        }

        //- The distance to the underlying searchableSurface
        scalar distance() const
        {
            return distance_;
        }

        //- The underlying surface
        const meshedSurface& surface() const
        {
            return *isoSurfacePtr_;
        }

        //- The underlying surface
        meshedSurface& surface()
        {
            return *isoSurfacePtr_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return labelList::null();
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            return isoSurfacePtr_->meshCells();
        }


    // Interpolate

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> interpolate
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


    // Output

        //- Print information
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
