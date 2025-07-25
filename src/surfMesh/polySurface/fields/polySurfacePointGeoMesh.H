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
    Foam::polySurfacePointGeoMesh

Description
    The polySurface GeoMesh for point fields.

\*---------------------------------------------------------------------------*/

#ifndef Foam_polySurfacePointGeoMesh_H
#define Foam_polySurfacePointGeoMesh_H

#include "polySurface.H"
#include "GeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class polySurfacePointGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class polySurfacePointGeoMesh
:
    public GeoMesh<polySurface>
{
public:

    // Constructors

        //- Construct from polySurface reference
        explicit polySurfacePointGeoMesh(const polySurface& mesh)
        :
            GeoMesh<polySurface>(mesh)
        {}


    // Static Functions

        //- The geometric (internal) size - number of points.
        //  Method name expected by GeoMesh interface
        static label size(const polySurface& mesh)
        {
            return mesh.nPoints();
        }

        //- The geometric boundary size - not used.
        //  Method name expected by GeoMesh interface
        static label boundary_size(const polySurface& mesh) noexcept
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

        // FUTURE?
        // //- Return the object registry for point fields
        // const objectRegistry& thisDb() const
        // {
        //     return mesh_.pointData();
        // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
