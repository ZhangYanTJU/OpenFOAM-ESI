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
    Check field availability for last time.
    Done to avoid mapping 'undefined' when a field only exists as time 0.

Requires
    readFields.H (for the checkData function)

\*---------------------------------------------------------------------------*/

// Initially all possible objects that are available at the final time
List<wordHashSet> availableRegionObjectNames(meshes.size());
List<HashTable<wordHashSet>> availableFaRegionObjectNames(meshes.size());

forAll(meshes, regioni)
{
    const auto& mesh = meshes[regioni];

    IOobjectList objects;
    HashTable<IOobjectList> faObjects;

    if (doConvertFields && !timeDirs.empty())
    {
        // const word& checkTimeDir = timeDirs.back().name();

        // List of volume mesh objects for this instance
        objects = IOobjectList(mesh, timeDirs.back().name());

        if (fieldSelector)
        {
            objects.filterObjects(fieldSelector);
        }

        // Remove "*_0" restart fields
        objects.prune_0();

        if (!doPointValues)
        {
            // Prune point fields if disabled
            objects.filterClasses(Foam::fieldTypes::is_point, true);
        }


        // The finite-area regions and their fields for this volume region
        // and instance
        faObjects.reserve(areaRegionNames.size());
        for (const word& areaName : areaRegionNames)
        {
            IOobjectList objs
            (
                faMesh::Registry(mesh),
                timeDirs.back().name(),
                polyMesh::regionName(areaName),
                IOobjectOption::NO_REGISTER
            );

            if (fieldSelector)
            {
                objs.filterObjects(fieldSelector);
            }

            // Remove "*_0" restart fields
            objs.prune_0();
            faObjects(areaName) = std::move(objs);
        }
    }

    // Volume fields
    if (!objects.empty())
    {
        wordList objectNames(objects.sortedNames());

        // Check availability for all times...
        checkData
        (
            mesh.thisDb(),
            timeDirs,
            objectNames
        );

        availableRegionObjectNames[regioni] = objectNames;
    }

    // Area fields
    for (const word& areaName : areaRegionNames)
    {
        wordList objectNames;

        if (const auto iter = faObjects.cfind(areaName); iter.good())
        {
            objectNames = iter.val().sortedNames();

            // Check availability for all times...
            checkData
            (
                faMesh::Registry(mesh),
                timeDirs,
                objectNames,
                polyMesh::regionName(areaName)
            );
        }

        availableFaRegionObjectNames[regioni].emplace_set
        (
            areaName,
            std::move(objectNames)
        );
    }
}


// ************************************************************************* //
