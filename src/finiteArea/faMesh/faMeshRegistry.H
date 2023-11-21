/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::faMeshRegistry

Description

    polyMesh
    |
    - finite-area
      - area1
      - area2

SourceFiles
    faMeshRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshRegistry_H
#define Foam_faMeshRegistry_H

#include "MeshObject.H"
#include "objectRegistry.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMeshRegistry
:
    public MeshObject<polyMesh, Foam::UpdateableMeshObject, faMeshRegistry>
{
    // Data Types

        //- The MeshObject type
        typedef
            MeshObject<polyMesh, Foam::UpdateableMeshObject, faMeshRegistry>
            meshObjectType;


    // Private Member Data

        //- The prefix to local: %finite-area
        static const word prefix_;

        //- The finite-area objects (meshes). Name: "finite-area"
        //  Anchored on the polyMesh parent
        objectRegistry objects_;


public:

    //- Runtime type information.
    TypeName("faMeshRegistry");

    //- Return the prefix to local: %finite-area
    static const word& prefix() noexcept;


    // Generated Methods

        //- No copy construct
        faMeshRegistry(const faMeshRegistry&) = delete;

        //- No copy assignment
        void operator=(const faMeshRegistry&) = delete;


    // Constructors

        //- Construct as singleton within polyMesh registry
        explicit faMeshRegistry(const polyMesh& mesh);


    // Database

        //- Return the object registry
        const objectRegistry& thisDb() const noexcept
        {
            Info<< "thisDb: " << objects_.name() << " = "
                << objects_.names() << nl;
            Info<< "parent: " << objects_.name() << " = "
                << objects_.parent().names() << nl;

            return objects_;
        }

        //- The polyMesh reference
        const polyMesh& mesh() const noexcept
        {
            return meshObjectType::mesh();
        }


    // Topological Change

        //- Is mesh moving
        bool moving() const { return meshObjectType::mesh().moving(); }

        //- Update after mesh motion
        virtual bool movePoints();

        //- Update after topo change
        virtual void updateMesh(const mapPolyMesh& mpm);


    // Write

        //- Write all objects
        virtual bool write(const bool writeOnProc = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
