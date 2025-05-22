/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021,2025 OpenCFD Ltd.
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
    Foam::surfPointGeoMesh

Description
    The surfMesh GeoMesh for point fields.

    Similar to surfGeoMesh, but refers to the surface points.

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfPointGeoMesh_H
#define Foam_surfPointGeoMesh_H

#include "GeoMesh.H"
#include "surfMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfPointGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class surfPointGeoMesh
:
    public GeoMesh<surfMesh>
{
public:

    // Constructors

        //- Construct from surfMesh reference
        explicit surfPointGeoMesh(const surfMesh& mesh)
        :
            GeoMesh<surfMesh>(mesh)
        {}


    // Static Functions

        //- The geometric (internal) size - number of points.
        //  Method name expected by GeoMesh interface
        static label size(const surfMesh& mesh)
        {
            return mesh.nPoints();
        }

        //- The geometric boundary size - not used.
        //  Method name expected by GeoMesh interface
        static label boundary_size(const surfMesh& mesh) noexcept
        {
            return 0;
        }


    // Member Functions

        //- The geometric (internal) size - number of points.
        //  Method name expected by GeoMesh interface
        label size() const
        {
            return size(mesh_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
