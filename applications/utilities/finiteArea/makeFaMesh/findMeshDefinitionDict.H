/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Search for the appropriate faMeshDefinition dictionary...

    Expects to find a faMeshDefinition file in a location specified
    by the command-line option, or one of the standard locations.

    For default area region (region0):
      (1) system/finite-area/faMeshDefinition
      (2) system/faMeshDefinition  [legacy location - v2312 and earlier]

    For a named area region:
      (1) system/finite-area/<area-name>/faMeshDefinition
      (2) system/finite-area/faMeshDefinition.<area-name>
          [symmetric with faOptions]

Required Classes
    - Foam::polyMesh
    - Foam::IOdictionary

Required Variables
    - regionName [word]  (the polyMesh region)
    - areaRegionName [word] (the areaMesh region)
    - args [argList]
    - runTime [Time]

Provided Variables
    - faMeshDefinitionPtr [autoPtr<IOdictionary>]

If the dictionary cannot be found, exits with an error message or
reports a warning (dry-run mode)

\*---------------------------------------------------------------------------*/

const word dictName("faMeshDefinition");

autoPtr<IOdictionary> faMeshDefinitionPtr;

{
    // Exit unless dry-run?
    const bool exitOnFailure = (!args.dryRun());

    const word& regionDir = polyMesh::regionName(regionName);
    const word& areaRegionDir = polyMesh::regionName(areaRegionName);

    // Canonical location
    fileName dictPath
    (
        runTime.system()/regionDir
      / faMesh::prefix()/areaRegionDir/dictName
    );


    // Dictionary specified on the command-line ...
    const bool specified = args.readIfPresent("dict", dictPath);

    if (specified)
    {
        if (isDir(dictPath))
        {
            dictPath /= dictName;
        }
    }


    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER,
        true    // is globalObject
    );

    refPtr<IOobject> meshDictIOPtr;

    bool foundIOobject = meshDictIO.typeHeaderOk<IOdictionary>(true);

    // For reporting any alternative locations
    dictPath.clear();

    if (!foundIOobject && !specified)
    {
        if (!areaRegionDir.empty())
        {
            // Alternative location
            // - system/finite-area/faMeshDefinition.<area-name>

            dictPath =
            (
                runTime.system()/regionDir
              / faMesh::prefix()/IOobject::groupName(dictName, areaRegionDir)
            );
        }
        else
        {
            // legacy location (v2312 and earlier)
            // - system/faMeshDefinition

            dictPath = runTime.system()/regionDir/dictName;
        }

        if (!dictPath.empty())
        {
            auto ioptr = autoPtr<IOobject>::New
            (
                dictPath,
                runTime,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER,
                true    // is globalObject
            );

            foundIOobject =
            (
                ioptr && ioptr->typeHeaderOk<IOdictionary>(true)
            );

            if (foundIOobject)
            {
                // Use if found
                meshDictIOPtr = std::move(ioptr);
            }

            // Generally retain dictPath for error messages,
            // but don't mention the really old legacy location
            if (areaRegionDir.empty())
            {
                dictPath.clear();
            }
        }
    }

    // Alternative location or regular location?
    if (!meshDictIOPtr)
    {
        meshDictIOPtr.cref(meshDictIO);
    }
    const auto& io = meshDictIOPtr();


    // Handle final success or failure
    if (foundIOobject)
    {
        Info<< "----------------" << nl
            << "Using area-mesh definition";
        if (!areaRegionDir.empty())
        {
            Info<< " [" << areaRegionName << "]";
        }

        Info<< nl
            << "    " << io.objectRelPath() << nl << endl;

        faMeshDefinitionPtr = autoPtr<IOdictionary>::New(io);
    }
    else
    {
        // Failure

        if (exitOnFailure)
        {
            auto& err =
                FatalErrorInFunction
                << "Did not find area-mesh definition";

            if (!areaRegionDir.empty())
            {
                err<< " [" << areaRegionName << "]";
            }

            err << nl
                << "    " << io.objectRelPath() << nl;

            if (!dictPath.empty())
            {
                err << "    " << dictPath << nl;
            }

            FatalError<< exit(FatalError);
        }
        else
        {
            // dry-run: warning
            auto& err =
                Warning << nl
                << "----------------" << nl
                << "Did not find area-mesh definition";

            if (!areaRegionDir.empty())
            {
                err << " [" << areaRegionName << "]";
            }

            err << nl
                << "    " << io.objectRelPath() << nl;

            if (!dictPath.empty())
            {
                err << "    " << dictPath << nl;
            }
        }
    }
}


// ************************************************************************* //
