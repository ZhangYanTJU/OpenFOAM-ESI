/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
------------------------------------------------------------------------------
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

#include "faOptions.H"
#include "faMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fa
    {
        defineTypeNameAndDebug(options, 0);
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Create IO object if dictionary is present
// - check finite-area locations
Foam::IOobject createIOobject
(
    const Foam::polyMesh& mesh,
    const Foam::word& baseName,  // eg, faOptions
    const Foam::word& areaName
)
{
    using namespace Foam;

    // eg, faOptions, faOptions.<area-name> etc
    const word lookupName
    (
        IOobject::groupName(baseName, polyMesh::regionName(areaName))
    );

    IOobject io
    (
        lookupName,
        mesh.time().constant(),
        // located under finite-area
        faMesh::Registry(mesh),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::REGISTER
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt(IOobjectOption::READ_MODIFIED);
    }
    else
    {
        // Check if faOptions, faOptions.<area-name> file is in system
        io.instance() = mesh.time().system();

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            io.readOpt(IOobjectOption::READ_MODIFIED);
        }
        else
        {
            io.readOpt(IOobjectOption::NO_READ);
        }
    }


    if (!io.isAnyRead() && polyMesh::regionName(areaName).empty())
    {
        // Check legacy location (default area region only)
        // - registered on polyMesh

        IOobject legacy
        (
            baseName,  // eg, faOptions
            mesh.time().constant(),
            mesh,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::REGISTER
        );

        if (legacy.typeHeaderOk<IOdictionary>(true))
        {
            legacy.readOpt(IOobjectOption::READ_MODIFIED);
        }
        else
        {
            // Check if the faOptions file is in system
            legacy.instance() = mesh.time().system();

            if (legacy.typeHeaderOk<IOdictionary>(true))
            {
                legacy.readOpt(IOobjectOption::READ_MODIFIED);
            }
            else
            {
                legacy.readOpt(IOobjectOption::NO_READ);
            }
        }

        if (legacy.isAnyRead())
        {
            Info<< "Creating finite-area options from "
                << legacy.instance()/legacy.name()
                << " (legacy location)" << nl << endl;

            return legacy;
        }
    }

    if (io.isAnyRead())
    {
        Info<< "Creating finite-area options from "
            << io.instance()/faMesh::prefix()/io.name() << nl << endl;
    }

    return io;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::options::options
(
    const fvMesh& mesh,
    const IOobject& io,
    const word& defaultAreaName
)
:
    IOdictionary(io),
    fa::optionList(mesh, *this, defaultAreaName)
{}


Foam::fa::options::options
(
    const fvMesh& mesh,
    const word& defaultAreaName
)
:
    IOdictionary(createIOobject(mesh, typeName, defaultAreaName)),
    fa::optionList(mesh, *this, defaultAreaName)
{}


Foam::fa::options& Foam::fa::options::New
(
    const fvMesh& mesh,
    const word& defaultAreaName
)
{
    // eg, faOptions, faOptions.<area-name> etc
    const word lookupName
    (
        IOobject::groupName(typeName, polyMesh::regionName(defaultAreaName))
    );

    // Registered under finite-area?
    auto* ptr = faMesh::Registry(mesh).getObjectPtr<fa::options>(lookupName);

    if (!ptr && polyMesh::regionName(defaultAreaName).empty())
    {
        // Legacy location?
        // Registered under polyMesh. region0 area only!
        ptr = mesh.thisDb().getObjectPtr<fa::options>(typeName);

        if (ptr)
        {
            InfoInFunction
                << "Retrieved  " << typeName
                << " from polyMesh " << mesh.name()
                << " (legacy location)" << endl;
        }
    }

    if (!ptr)
    {
        DebugInFunction
            << "Constructing " << lookupName
            << " for region " << mesh.name()
            << " : " << polyMesh::regionName(defaultAreaName) << endl;

        ptr = new fa::options
        (
            mesh,
            createIOobject(mesh, typeName, defaultAreaName),
            defaultAreaName
        );

        regIOobject::store(ptr);
    }

    return *ptr;
}


bool Foam::fa::options::read()
{
    if (IOdictionary::regIOobject::read())
    {
        fa::optionList::read(*this);
        return true;
    }

    return false;
}


// ************************************************************************* //
