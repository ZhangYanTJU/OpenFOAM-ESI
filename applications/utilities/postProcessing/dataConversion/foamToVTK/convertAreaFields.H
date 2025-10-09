/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for converting finite-area - included by foamToVTK.

    Typedef vtkWriterType_areaMesh is declared in writeAreaFields.H

\*---------------------------------------------------------------------------*/

//
// Finite-area mesh and fields - need not exist
//

// No subsetting!
if (doFiniteArea && !areaRegionNames.empty())
{
    using reportFields = foamToVtkReportFields;

    // The (region) polyMesh being used. No subsetting possible
    const auto& basePolyMesh = meshProxy.baseMesh();

    for (const word& areaName : areaRegionNames)
    {
        const bool isDefaultRegion(polyMesh::regionName(areaName).empty());

        // CAUTION
        // If we want to have constant access to the HashTable:
        //
        // (SEGFAULT)
        // const auto& faObjs = faObjects.lookup(areaName, IOobjectList());
        //
        // Use an empty fallback to avoid binding to a temporary:
        //
        // const IOobjectList emptyObjectList;
        // const auto& faObjs = faObjects.lookup(areaName, emptyObjectList);

        // Since we do not need the area fields afterwards,
        // just move them out from the HashTable

        IOobjectList faObjs;
        if (auto iter = faObjects.find(areaName); iter.good())
        {
            faObjs = std::move(iter.val());
        }

        const label nAreaFields = faObjs.count(Foam::fieldTypes::is_area);

        autoPtr<faMesh> faMeshPtr;

        if (nAreaFields || withMeshIds)
        {
            faMeshPtr = faMesh::TryNew(areaName, basePolyMesh);
        }

        if (!faMeshPtr)
        {
            if (!isDefaultRegion)
            {
                // Report any area region specified but missing
                // - silently ignore region0
                Info<< "No area-mesh [" << polyMesh::regionName(areaName)
                    << "] on volume-region ["
                    << basePolyMesh.regionName() << "]" << endl;
            }
            continue;
        }
        const auto& areaMesh = faMeshPtr();

        Info<< "Using area-mesh [" << polyMesh::regionName(areaName)
            << "] on volume-region ["
            << basePolyMesh.regionName() << "]" << endl;


        reportFields::area(Info, faObjs);

        const auto& pp = faMeshPtr->patch();

        vtkWriterType_areaMesh writer
        (
            pp,
            writeOpts,
            (
                outputDir/regionDir/"finite-area"
              / (
                   isDefaultRegion
                 ? fileName("finiteArea")
                 : fileName(areaName/areaName)
                )
              + timeDesc
            ),
            UPstream::parRun()
        );
        Info<< "    Area      : "
            << args.relativePath(writer.output()) << nl;

        writer.beginFile(areaMesh.name());

        writer.writeTimeValue(timeValue);
        writer.writeGeometry();

        // Optionally with (cellID, patchID, faceLabels, procID) fields
        writer.beginCellData
        (
            (withMeshIds ? 3 + (writer.parallel() ? 1 : 0) : 0)
          + nAreaFields
        );

        if (withMeshIds)
        {
            const globalIndex procAddr(areaMesh.nFaces());

            // Use global indexed values for the 'cell' ids
            writer.writeCellData("cellID", identity(procAddr.range()));

            // The patch ids can also be quite useful
            const polyBoundaryMesh& pbm = areaMesh.mesh().boundaryMesh();

            labelList patchIds
            (
                pbm.patchID(areaMesh.faceLabels())
            );

            writer.writeCellData("patchID", patchIds);

            // Use proc-local data for faceLabels
            // (confusing enough already without renumbering)
            writer.writeCellData("faceLabels", areaMesh.faceLabels());

            writer.writeProcIDs();  // parallel only
        }

        writeAllAreaFields
        (
            writer,
            areaMesh,
            faObjs,
            true // syncPar
        );

        fileName outputName(writer.output());

        writer.close();

        if (UPstream::master())
        {
            // Add to file-series and emit as JSON

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }
    }
}


// ************************************************************************* //
