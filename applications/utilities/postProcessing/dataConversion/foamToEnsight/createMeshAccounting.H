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
    Additional mesh accounting (Ensight)

\*---------------------------------------------------------------------------*/

// Cases and meshes per volume region
PtrList<ensightCase> ensightCases(regionNames.size());
PtrList<ensightMesh> ensightMeshes(regionNames.size());

// Per volume region can have multiple finite-area regions
List<PtrDynList<faMesh>> meshesFa(regionNames.size());
List<PtrDynList<ensightCase>> ensightCasesFa(regionNames.size());
List<PtrDynList<ensightFaMesh>> ensightMeshesFa(regionNames.size());

{
    forAll(regionNames, regioni)
    {
        const fvMesh& mesh = meshes[regioni];

        const word& regionName = regionNames[regioni];
        const word& regionDir = polyMesh::regionName(regionName);

        fileName ensCasePath(outputDir);
        word ensCaseName(args.globalCaseName());

        if (!regionDir.empty())
        {
            ensCaseName = regionName;
            ensCasePath /= regionName;

            // Handle very rare naming collision with Ensight directories
            if (regionName == "data")
            {
                ensCasePath += "-region";
            }
        }

        // Volume mesh
        ensightMeshes.set
        (
            regioni,
            new ensightMesh(mesh, writeOpts)
        );
        ensightMeshes[regioni].verbose(optVerbose);

        // New ensight case file, initialize header etc.
        ensightCases.set
        (
            regioni,
            new ensightCase(ensCasePath, ensCaseName, caseOpts)
        );

        if (!doFiniteArea)
        {
            continue;
        }

        // Note: not every volume region is guaranteed to have
        // any or all area region(s)

        meshesFa[regioni].reserve_exact(areaRegionNames.size());

        for (const word& areaName : areaRegionNames)
        {
            auto faMeshPtr = faMesh::TryNew(areaName, mesh);

            if (faMeshPtr)
            {
                meshesFa[regioni].push_back(std::move(faMeshPtr));
            }
        }

        // Setup the ensight components
        const label nAreas = meshesFa[regioni].size();
        ensightCasesFa[regioni].reserve_exact(nAreas);
        ensightMeshesFa[regioni].reserve_exact(nAreas);

        for (const faMesh& areaMesh : meshesFa[regioni])
        {
            // Use regionName() - automatically filters for defaultRegion
            const word& areaName = areaMesh.regionName();

            if (areaName.empty())
            {
                // single-region
                ensightCasesFa[regioni].emplace_back
                (
                    ensCasePath/"finite-area",
                    "finite-area",
                    caseOpts
                );
            }
            else
            {
                // multi-region
                ensightCasesFa[regioni].emplace_back
                (
                    ensCasePath/"finite-area"/areaName,
                    areaName,
                    caseOpts
                );
            }

            auto& ensFaMesh = ensightMeshesFa[regioni].emplace_back(areaMesh);
            ensFaMesh.verbose(optVerbose);
        }
    }
}


// ************************************************************************* //
