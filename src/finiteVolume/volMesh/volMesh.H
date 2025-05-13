/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::volMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

\*---------------------------------------------------------------------------*/

#ifndef Foam_volMesh_H
#define Foam_volMesh_H

#include "GeoMesh.H"
#include "fvMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class volMesh Declaration
\*---------------------------------------------------------------------------*/

class volMesh
:
    public GeoMesh<fvMesh>
{
public:

    // Constructors

        //- Construct from fvMesh reference
        explicit volMesh(const fvMesh& mesh)
        :
            GeoMesh<fvMesh>(mesh)
        {}


    // Static Functions

        //- The geometric (internal) mesh size - number of cells.
        //  Method name expected by GeoMesh interface
        static label size(const polyMesh& mesh) noexcept
        {
            return mesh.nCells();
        }

        //- The geometric boundary size - number of boundary faces.
        //  Method name expected by GeoMesh interface
        static label boundary_size(const polyMesh& mesh) noexcept
        {
            return mesh.nBoundaryFaces();
        }


    // Member Functions

        //- The geometric (internal) mesh size - number of cells.
        //  Method name expected by GeoMesh interface
        label size() const noexcept
        {
            return size(mesh_);
        }

        //- Field of cell centres
        const volVectorField& C() const
        {
            return mesh_.C();
        }
};


/*---------------------------------------------------------------------------*\
                          Class isVolMesh Declaration
\*---------------------------------------------------------------------------*/

//- Template argument for type volMesh.
template<class> struct isVolMesh : std::false_type {};

template<> struct isVolMesh<volMesh> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
