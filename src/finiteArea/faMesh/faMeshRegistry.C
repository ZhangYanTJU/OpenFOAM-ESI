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

\*---------------------------------------------------------------------------*/

#include "faMeshRegistry.H"
#include "faMesh.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMeshRegistry, 0);
}

const Foam::word Foam::faMeshRegistry::prefix_("finite-area");

const Foam::word& Foam::faMeshRegistry::prefix() noexcept
{
    return prefix_;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

#if 0
const Foam::faMeshRegistry*
Foam::faMeshRegistry::hasRegistry(const polyMesh& mesh)
{
    return mesh.objectRegistry::cfindObject<faMeshRegistry>
    (
        faMeshRegistry::typeName
    );
}


const Foam::faMeshRegistry&
Foam::faMeshRegistry::registry(const polyMesh& mesh)
{
    // Force create on access
    return faMeshRegistry::New(mesh);
}


Foam::UPtrList<const Foam::faMesh>
Foam::faMeshRegistry::meshes(const polyMesh& mesh)
{
    const auto* subreg =
        mesh.objectRegistry::cfindObject<faMeshRegistry>
        (
            faMeshRegistry::typeName
        );

    if (subreg)
    {
        return subreg->csorted<faMesh>();
    }

    return UPtrList<const faMesh>();
}


Foam::label Foam::faMeshRegistry::numMeshes(const polyMesh& mesh)
{
    const auto* subreg =
        mesh.objectRegistry::cfindObject<faMeshRegistry>
        (
            faMeshRegistry::typeName
        );

    if (subreg)
    {
        return subreg->count<faMesh>();
    }

    return 0;
}

#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshRegistry::faMeshRegistry(const polyMesh& mesh)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMeshRegistry>(mesh),
    objects_
    (
        IOobject
        (
            faMeshRegistry::prefix(),
            mesh.time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        )
    )
{
    Info<< "time: " << mesh.time().sortedNames() << nl;
    Info<< "mesh: " << mesh.sortedNames() << nl;

    Info<< "Construct registry..." << endl;

    {
        HashTable<const regIOobject*> objs =
            mesh.lookupClass<const regIOobject>();

        forAllConstIters(objs, iter)
        {
            Info<< "    " << iter.key() << " -> " << iter.val()->type() << nl;
        }

        Info<< "mesh: " << mesh.sortedNames() << nl;
    }

    {
        HashTable<const regIOobject*> objs =
            objects_.lookupClass<const regIOobject>();

        forAllConstIters(objs, iter)
        {
            Info<< "    " << iter.key() << " -> " << iter.val()->type() << nl;
        }

        Info<< "registry: " << objects_.sortedNames() << nl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// const Foam::polyMesh Foam::faMeshRegistry::mesh() const
// {
//     return Mesh
//     MeshObject<polyMesh, Foam::UpdateableMeshObject, faMeshRegistry>(mesh),
//     objects_
//     (
//         IOobject
//         (
//             faMeshRegistry::prefix(),
//             mesh.time().timeName(),
//             mesh.db(),
//             IOobjectOption::NO_READ,
//             IOobjectOption::AUTO_WRITE,
//             IOobjectOption::REGISTER
//         )
//     )

bool Foam::faMeshRegistry::movePoints()
{
    UPtrList<faMesh> list = objects_.sorted<faMesh>();

    for (faMesh& m : list)
    {
        m.movePoints();
    }

    return true;
}


void Foam::faMeshRegistry::updateMesh(const mapPolyMesh& mpm)
{
    UPtrList<faMesh> list = objects_.sorted<faMesh>();

    for (faMesh& m : list)
    {
        m.updateMesh(mpm);
    }
}


bool Foam::faMeshRegistry::write(const bool writeOnProc) const
{
    UPtrList<const faMesh> list = objects_.csorted<faMesh>();

    for (const faMesh& m : list)
    {
        m.write(writeOnProc);
    }

    return true;
}


// ************************************************************************* //
