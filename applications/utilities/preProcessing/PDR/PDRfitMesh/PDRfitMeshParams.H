/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PDRfitMeshParams

Description
    Parameters for PDRfitMesh

SourceFiles
    PDRfitMeshParams.C

\*---------------------------------------------------------------------------*/

#ifndef PDRfitMeshParams_H
#define PDRfitMeshParams_H

#include "labelList.H"
#include "scalarList.H"
#include "dictionary.H"
#include "optionalData.H"
#include "MinMax.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PDRfitMeshParams Declaration
\*---------------------------------------------------------------------------*/

class PDRfitMeshParams
{
    // Private Member Functions

        //- Read optional bounds
        void readBounds(const dictionary& dict);


public:

    // Data Members

        //- Optional user bounds (xmin, xmax, ymin, ymax, zmin, zmax)
        Vector<optionalData<scalar>> minBounds, maxBounds;

        //- Ground datum (hard zmin)
        scalar ground = -1e16;


    // Finding planes

        //- Ignore faces smaller than this multiple of cell face area
        scalar minFaceArea  = 5.0;

        //- Only fit a plane if the face area at the coordinate is at
        //- least this times the cell face area
        scalar minAreaRatio = 20.0;

        //- Size of search zones for face areas will be this * the
        //- cell width.
        //  Faces closer than this zone size may be grouped together
        scalar areaWidthFactor = 0.7;

        // Very long zones produce bad outer boundary shape from
        // makePDRmeshBlocks, so we subdivide a zone if its length is
        // greater than this * the height of the area of cuboid vells
        scalar maxZoneToHeight = 2.0;


    // Choosing cell width

        // Minimum number of cells in any direction
        label nCellsMin = 10;

        // Width of obstacles must be less than this * cell width to
        // added into subgrid length
        scalar widthFactor = 1.0;

        // Optimal average number of obstacles per cell
        scalar obsPerCell = 2.0;

        //- Maximum cellWidth when auto-sizing
        scalar maxCellWidth = 1e16;

        //- Do not use a cellWidth more than this times the initial estimate
        scalar maxWidthEstimate = 5.0;

        //- Converged if optimised width changes by less than this amount
        scalar maxWidthRatio = 1.2;

        //- Maximum iterations in optimising cellWidth
        label maxIterations = 5;


    // Outer region

        //- Fraction of rectangular cell layers on each side of the
        //- central region
        scalar nEdgeLayers = 5;

        //- Cell size ratio - should be same as used for PDRblockMesh
        scalar outerRatio = 1.2;

//FUTURE?  // Ratio of the outer extension distance to the average
//FUTURE?  // radius of the congested region
//FUTURE?  scalar outerRadius = 20.0;


    // Constructors

        //- Default construct
        PDRfitMeshParams() = default;

        //- Construct and read dictionary
        explicit PDRfitMeshParams(const dictionary& dict);


    // Member Functions

        //- Read program parameters from dictionary
        void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
