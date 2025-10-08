/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

Application
    makeFaMesh

Description
    Check a finite-area mesh

Original Authors
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "faMesh.H"
#include "faMeshTools.H"
#include "polyMesh.H"
#include "areaFaMesh.H"
#include "edgeFaMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "processorFaPatch.H"
#include "foamVtkIndPatchWriter.H"
#include "foamVtkLineWriter.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check a finite-area mesh"
    );

    argList::addBoolOption
    (
        "write-vtk",
        "Write mesh as a vtp (vtk) file for display or debugging"
    );

    argList::addOption
    (
        "geometryOrder",
        "N",
        "Test different geometry order - experimental!!",
        true  // Advanced option
    );

    #include "addRegionOption.H"
    #include "addAllFaRegionOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Handle -all-area-regions, -area-regions, -area-region
    #include "getAllFaRegionOptions.H"

    // ------------------------------------------------------------------------

    #include "createNamedPolyMesh.H"

    int geometryOrder(1);
    if (args.readIfPresent("geometryOrder", geometryOrder))
    {
        Info<< "Setting faMesh::geometryOrder = " << geometryOrder << nl
            << "(experimental)" << nl << endl;

        faMesh::geometryOrder(geometryOrder);
    }

    for (const word& areaName : areaRegionNames)
    {
        Info<< "Create faMesh";
        if (!polyMesh::regionName(areaName).empty())
        {
            Info<< " [" << areaName << "]";
        }
        Info<< " for time = " << runTime.timeName() << nl;

        autoPtr<faMesh> faMeshPtr(faMesh::TryNew(areaName, mesh));

        if (!faMeshPtr)
        {
            Info<< "    ...failed to create area-mesh";
            if (!polyMesh::regionName(areaName).empty())
            {
                Info<< " [" << areaName << "]";
            }
            Info<< endl;
            continue;
        }
        else
        {
            Info<< endl;
        }

        const auto& aMesh = faMeshPtr();

        // Mesh information (verbose)
        faMeshTools::printMeshChecks(aMesh);

        if (args.found("write-vtk"))
        {
            #include "faMeshWriteVTK.H"
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
