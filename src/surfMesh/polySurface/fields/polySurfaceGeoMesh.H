/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::polySurfaceGeoMesh

Description
    The polySurface GeoMesh for face fields.

\*---------------------------------------------------------------------------*/

#ifndef Foam_polySurfaceGeoMesh_H
#define Foam_polySurfaceGeoMesh_H

#include "polySurface.H"
#include "GeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polySurfaceGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class polySurfaceGeoMesh
:
    public GeoMesh<polySurface>
{
public:

    // Constructors

        //- Construct from polySurface reference
        explicit polySurfaceGeoMesh(const polySurface& mesh)
        :
            GeoMesh<polySurface>(mesh)
        {}


    // Static Functions

        //- The geometric (internal) size - number of faces.
        //  Method name expected by GeoMesh interface
        static label size(const polySurface& mesh)
        {
            return mesh.nFaces();
        }

        //- The geometric boundary size - not used.
        //  Method name expected by GeoMesh interface
        static label boundary_size(const polySurface& mesh) noexcept
        {
            return 0;
        }


    // Member Functions

        //- The geometric (internal) size - number of faces.
        //  Method name expected by GeoMesh interface
        label size() const
        {
            return size(mesh_);
        }

        // FUTURE?
        // //- Return the object registry for face fields
        // const objectRegistry& thisDb() const
        // {
        //     return mesh_.faceData();
        // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
