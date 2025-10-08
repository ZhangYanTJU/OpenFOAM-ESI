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
    VTK output of faMesh with some geometric or debug fields

\*---------------------------------------------------------------------------*/

// The base name for output files
word vtkBaseFileName(aMesh.regionName());
if (polyMesh::regionName(vtkBaseFileName).empty())
{
    vtkBaseFileName = "finiteArea";
}

Info<< nl;
Info<< "Write faMesh in vtk format:" << nl;

{
    // finiteArea - faces
    vtk::uindirectPatchWriter writer
    (
        aMesh.patch(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.time().globalPath() / vtkBaseFileName
        )
    );

    writer.writeGeometry();

    globalIndex procAddr(aMesh.nFaces());
    labelList cellIDs;

    if (UPstream::master())
    {
        cellIDs.resize(procAddr.totalSize());
        for (const labelRange& range : procAddr.ranges())
        {
            auto slice = cellIDs.slice(range);
            slice = identity(range);
        }
    }

    // CellData
    writer.beginCellData(4);
    writer.writeProcIDs();
    writer.write("cellID", cellIDs);
    writer.write("area", aMesh.S().field());
    writer.write("normal", aMesh.faceAreaNormals());

    // PointData
    writer.beginPointData(1);
    writer.write("normal", aMesh.pointAreaNormals());

    Info<< "    " << writer.output().name() << nl;
}

{
    // finiteArea - edges
    vtk::lineWriter writer
    (
        aMesh.points(),
        aMesh.edges(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.time().globalPath() / (vtkBaseFileName + "-edges")
        )
    );

    writer.writeGeometry();

    // CellData
    writer.beginCellData(4);
    writer.writeProcIDs();
    {
        Field<scalar> fld
        (
            // Use primitive patch order
            faMeshTools::flattenEdgeField(aMesh.magLe(), true)
        );
        writer.write("magLe", fld);
    }

    // PointData
    writer.beginPointData(1);
    writer.write("normal", aMesh.pointAreaNormals());

    Info<< "    " << writer.output().name() << nl;
}

{
    const Field<vector> edgeCentres
    (
        // Use primitive patch order
        faMeshTools::flattenEdgeField(aMesh.edgeCentres(), true)
    );

    // finiteArea - edgeCentres
    // (no other convenient way to display vectors on the edges)
    vtk::lineWriter writer
    (
        edgeCentres,
        edgeList::null(),
        // vtk::formatType::INLINE_ASCII,
        fileName
        (
            aMesh.time().globalPath() / (vtkBaseFileName + "-edgesCentres")
        )
    );

    writer.writeGeometry();

    // PointData
    writer.beginPointData(3);
    writer.writeProcIDs();  // Unfortunately cannot threshold on points
    {
        Field<vector> fld
        (
            // Use primitive patch order
            faMeshTools::flattenEdgeField(aMesh.Le(), true)
        );
        writer.write("Le", fld);
    }
    {
        Field<vector> fld
        (
            // Use primitive patch order
            faMeshTools::flattenEdgeField(aMesh.edgeAreaNormals(), true)
        );
        writer.write("normal", fld);
    }

    Info<< "    " << writer.output().name() << nl;
}


// ************************************************************************* //
