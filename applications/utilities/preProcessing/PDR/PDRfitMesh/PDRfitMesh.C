/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

Applications
    PDRfitMesh

Description
    Scans extents of obstacles to establish reasonable estimates
    for generating a PDRblockMeshDict.

SourceFiles
    PDRfitMesh.C

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Fstream.H"
#include "IOdictionary.H"

#include "PDRfitMeshScans.H"
#include "PDRsetFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Processes a set of geometrical obstructions to determine"
        " reasonable estimates for generating a PDRblockMeshDict"
    );
    argList::noParallel();
    argList::noFunctionObjects();

    argList::addOption("dict", "file", "Alternative PDRsetFieldsDict");
    argList::addOption("params", "file", "Alternative PDRfitMeshDict");
    argList::addBoolOption
    (
        "overwrite",
        "Overwrite existing system/PDRblockMeshDict"
    );

    argList::addBoolOption("verbose", "Increase verbosity");

    argList::addBoolOption
    (
        "dry-run",
        "Equivalent to -print-dict"
    );

    argList::addBoolOption
    (
        "print-dict",
        "Print PDRblockMeshDict equivalent and exit"
    );

    argList::addBoolOption
    (
        "write-vtk",
        "Write obstacles as VTK files"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool dryrun = args.found("dry-run");
    const bool printDict = args.found("print-dict");

    const word dictName("PDRsetFieldsDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary setFieldsDict(dictIO);

    const fileName& casepath = runTime.globalPath();

    // Program parameters (globals)
    pars.read(setFieldsDict);

    // Params for fitMesh
    // - like setSystemRunTimeDictionaryIO

    IOobject paramIO = IOobject::selectIO
    (
        IOobject
        (
            "PDRfitMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        args.getOrDefault<fileName>("params", "")
    );

    if (paramIO.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Using PDRfitMesh parameters from "
            << runTime.relativePath(paramIO.objectPath()) << nl
            << endl;
    }
    else
    {
        paramIO.readOpt() = IOobject::NO_READ;
        Info<< "No PDRfitMeshDict found, using defaults" << nl
            << endl;
    }


    IOdictionary fitParamsDict(paramIO, dictionary());

    PDRfitMeshParams fitParams(fitParamsDict);


    // Essential parameters

    scalar cellWidth = 0;

    if
    (
        !setFieldsDict.readIfPresent("cellWidth", cellWidth)
     && !fitParamsDict.readIfPresent("cellWidth", cellWidth)
    )
    {
        FatalErrorInFunction
            << "No cellWidth specified in any dictionary" << nl
            << exit(FatalError);
    }


    // Storage for obstacles and cylinder-like obstacles
    DynamicList<PDRobstacle> obstacles, cylinders;

    // Read in obstacles
    if (pars.legacyObsSpec)
    {
        PDRobstacle::legacyReadFiles
        (
            pars.obsfile_dir, pars.obsfile_names,
            boundBox::invertedBox,  // ie, no bounds
            obstacles,
            cylinders
        );
    }
    else
    {
        PDRobstacle::readFiles
        (
            pars.obsfile_dir, pars.obsfile_names,
            boundBox::invertedBox,  // ie, no bounds
            obstacles,
            cylinders
        );
    }


    //
    // Output names
    //

    IOobject outputIO
    (
        "PDRblockMeshDict",
        runTime.system(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false  // unregistered
    );


    enum class writeType { DRY_RUN, OKAY, OVERWRITE, NO_CLOBBER };

    writeType writable = writeType::OKAY;

    if (dryrun || printDict)
    {
        writable = writeType::DRY_RUN;
    }
    else if (isFile(outputIO.objectPath()))
    {
        if (args.found("overwrite"))
        {
            writable = writeType::OVERWRITE;
        }
        else
        {
            writable = writeType::NO_CLOBBER;
            InfoErr
                << nl
                << "File exists: "
                << runTime.relativePath(outputIO.objectPath()) << nl
                << "Move out of the way or specify -overwrite" << nl << nl
                << "Exiting" << nl << nl;

            return 1;
        }
    }
    else
    {
        writable = writeType::OKAY;
    }


    if (args.found("write-vtk"))
    {
        PDRobstacle::generateVtk(casepath/"VTK", obstacles, cylinders);
    }


    //
    // The fitting routines
    //

    // Collapse into a single list of obstacles
    obstacles.append(std::move(cylinders));

    PDRfitMeshScan::verbose(args.found("verbose"));

    Vector<PDRblock::gridControl> griding =
        PDRfitMeshScans().calcGriding(obstacles, fitParams, cellWidth);


    //
    // Output
    //

    if (writable == writeType::DRY_RUN)
    {
        InfoErr << nl;
        if (!printDict)
        {
            InfoErr
                << "dry-run: ";
        }
        InfoErr
            << "Displaying equivalent PDRblockMeshDict" << nl
            << nl;
    }
    else if (writable == writeType::OKAY)
    {
        InfoErr
            << nl
            << "Write "
            << runTime.relativePath(outputIO.objectPath()) << nl;
    }
    else if (writable == writeType::OVERWRITE)
    {
        InfoErr
            << nl
            << "Overwrite existing "
            << runTime.relativePath(outputIO.objectPath()) << nl;
    }
    // NO_CLOBBER already handled


    {
        autoPtr<Ostream> outputFilePtr;

        if
        (
            writable == writeType::OKAY
         || writable == writeType::OVERWRITE
        )
        {
            outputFilePtr.reset(new OFstream(outputIO.objectPath()));
        }

        Ostream& os = bool(outputFilePtr) ? *outputFilePtr : Info();

        outputIO.writeHeader(os, IOdictionary::typeName_());

        os.writeEntry
        (
            "expansion",
            PDRfitMeshScan::expansionName()
        );
        os  << nl;

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            griding[cmpt].writeDict(os, cmpt);
        }

        const dictionary* outerDict;

        if
        (
            (outerDict = setFieldsDict.findDict("outer")) != nullptr
         || (outerDict = fitParamsDict.findDict("outer")) != nullptr
        )
        {
            outerDict->writeEntry(os);
        }
        else
        {
            // Define our own "outer"
            os.beginBlock("outer");
            os.writeEntry("type", "none");
            os.endBlock();
        }

        IOobject::writeEndDivider(os);
    }

    Info<< nl << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
