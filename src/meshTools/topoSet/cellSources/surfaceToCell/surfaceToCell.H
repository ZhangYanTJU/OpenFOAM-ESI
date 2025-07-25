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
    Foam::surfaceToCell

Description
    A \c topoSetCellSource to select cells based on
    relation to a surface given by an external file.

    \c surfaceToCell can select:
    - all cells inside/outside/cut by a surface
    - all cells inside/outside surface (\c useSurfaceOrientation,
    requires closed surface)
    - cells with centre nearer than XXX to surface
    - cells with centre nearer than XXX to surface \b and with normal
      at nearest point to centre and cell-corners differing by
      more than YYY (i.e. point of high curvature)

    Operands:
    \table
      Operand    | Type       | Location
      input      | triSurface | $FOAM_CASE/constant/triSurface/\<file\>
      output     | cellSet    | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name            <name>;
        type            cellSet;
        action          <action>;

        // Mandatory entries
        source          surfaceToCell;
        file            <surfaceFileName>;
        outsidePoints
        (
            (<p1x> <p1y> <p1z>)
            (<p2x> <p2y> <p2z>)
            ...
        );
        includeCut      false;
        includeInside   false;
        includeOutside  true;
        nearDistance    0.5;
        curvature       1.0;

        // Optional entries
        useSurfaceOrientation false;
        fileType        stl;
        scale           2.0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: surfaceToCell          | word |  yes  | -
      file       | The surface "filename"              | word |  yes  | -
      outsidePoints | List of points that define outside of the surface <!--
                                                   --> | vectorList | yes | -
      includeCut | Flag to include cut cells           | bool |  yes  | -
      includeInside | Flag to include inside cells     | bool |  yes  | -
      includeOutside | Flag to include outside cells   | bool |  yes  | -
      useSurfaceOrientation | Flag to use inherently the orientation of <!--
                            --> the surface            | bool |  no   | false
      nearDistance | Near distance to the surface      | scalar | yes | -
      curvature    | Surface curvature                 | scalar | yes | -
      fileType     | The format of the surface file    | word   | no  | ""
      scale        | Surface scaling factor            | scalar | no  | -1
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    surfaceToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceToCell_H
#define Foam_surfaceToCell_H

#include "topoSetCellSource.H"
#include "refPtr.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class triSurface;
class triSurfaceSearch;

/*---------------------------------------------------------------------------*\
                           Class surfaceToCell Declaration
\*---------------------------------------------------------------------------*/

class surfaceToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of surface file
        const fileName surfName_;

        //- Points which are outside
        const pointField outsidePoints_;

        //- Include cut cells
        const bool includeCut_;

        //- Include inside cells
        const bool includeInside_;

        //- Include outside cells
        const bool includeOutside_;

        //- Determine inside/outside purely using geometric test
        //  (does not allow includeCut)
        const bool useSurfaceOrientation_;

        //- If > 0 : include cells with distance from cellCentre to surface
        //  less than nearDist.
        const scalar nearDist_;

        //- If > -1 : include cells with normals at nearest surface points
        //  varying more than curvature_.
        const scalar curvature_;

        //- triSurface to search on (can be external reference)
        const refPtr<triSurface> surfPtr_;

        //- Search engine on surface.
        const refPtr<triSurfaceSearch> querySurfPtr_;


    // Private Member Functions

        //- Find index of nearest triangle to point. Returns triangle or -1 if
        //  not found within search span.
        //  Cache result under pointi.
        static label getNearest
        (
            const triSurfaceSearch& querySurf,
            const label pointi,
            const point& pt,
            const vector& searchSpan,
            Map<label>& cache
        );

        //- Return true if surface normal of nearest points to vertices on
        //  cell differ from that on cell centre. Points cached in
        //  pointToNearest.
        bool differingPointNormals
        (
            const triSurfaceSearch& querySurf,
            const vector& span,
            const label celli,
            const label cellTriI,
            Map<label>& pointToNearest
        ) const;


        //- Depending on surface add to or delete from cellSet.
        void combine(topoSet& set, const bool add) const;

        //- Check values at construction time.
        void checkSettings() const;

        const triSurfaceSearch& querySurf() const
        {
            return *querySurfPtr_;
        }


public:

    //- Runtime type information
    TypeName("surfaceToCell");

    // Constructors

        //- Construct from components
        surfaceToCell
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const pointField& outsidePoints,
            const bool includeCut,
            const bool includeInside,
            const bool includeOutside,
            const bool useSurfaceOrientation,
            const scalar nearDist,
            const scalar curvature
        );

        //- Construct from components (supplied surface, surfaceSearch)
        surfaceToCell
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const triSurface& surf,
            const triSurfaceSearch& querySurf,
            const pointField& outsidePoints,
            const bool includeCut,
            const bool includeInside,
            const bool includeOutside,
            const bool useSurfaceOrientation,
            const scalar nearDist,
            const scalar curvature
        );

        //- Construct from dictionary
        surfaceToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        surfaceToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~surfaceToCell();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
