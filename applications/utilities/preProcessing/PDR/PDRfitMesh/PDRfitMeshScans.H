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
    Foam::PDRfitMeshScans

Description
    Scanning of obstacles in a multiple directions for PDRfitMesh

SourceFiles
    PDRfitMeshScans.C

\*---------------------------------------------------------------------------*/

#ifndef PDRfitMeshScans_H
#define PDRfitMeshScans_H

#include "PDRfitMeshScan.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class PDRobstacle;

/*---------------------------------------------------------------------------*\
                       Class PDRfitMeshScans Declaration
\*---------------------------------------------------------------------------*/

class PDRfitMeshScans
:
    public Vector<PDRfitMeshScan>
{
    // Private Data

        //- The evaluated sub-grid length
        scalar subgridLen_ = 0;


    // Private Member Functions

        //- Prepare limits, scan areas
        void prepare
        (
            const UList<PDRobstacle>& obstacles,
            const PDRfitMeshParams& fitParams,
            scalar cellWidth
        );

public:

    // Default Constructors


    // Member Functions

        //- Verbosity
        static inline bool verbose() noexcept
        {
            return PDRfitMeshScan::verbose();
        }

        //- Calculate the grid controls for the given obstacles
        //- and parameters
        Vector<PDRblock::gridControl>
        calcGriding
        (
            const UList<PDRobstacle>& obstacles,
            const PDRfitMeshParams& fitParams,
            scalar cellWidth
        );


    // Low-level functions

        //- The volume of the limits
        scalar volume() const;

        //- The volume of single cell of the respective step-width
        scalar cellVolume() const;

        //- Reset to initial state
        void reset();

        //- Define number of steps (and step-size) according to the
        //- current limits and the specified parameters
        void resize(const scalar cellWidth, const PDRfitMeshParams& pars);

        //- Adjust directional limits to accommodate obstacles
        void scanLimits(const UList<PDRobstacle>& obstacles);

        //- Populate with areas/positions of the obstacles,
        //- and check for sub-grid obstacles
        //  Sets the internal
        void scanAreas
        (
            const UList<PDRobstacle>& obstacles,
            const PDRfitMeshParams& fitParams,
            scalar cellWidth
        );

        //- Populate with areas/positions of the obstacles,
        //- and check for sub-grid obstacles
        //
        //  \return total subgrid lengths
        scalar scanAreas
        (
            const UList<PDRobstacle>& obstacles,
            const scalar minSubgridLen
        );

        //- Populate with areas/positions of the obstacles,
        //- without checks for sub-grid obstacles
        void scanAreas(const UList<PDRobstacle>& obstacles);

        //- Print information
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
