/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

#include "loadOrCreateMesh.H"
#include "faMesh.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "decomposedBlockData.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// Trimmed-down version of lookupAndCacheProcessorsPath
// with Foam::exists() check. No caching.

// Check for two conditions:
// - file has to exist
// - if collated the entry has to exist inside the file

// Note: bypass fileOperation::filePath(IOobject&) since has problems
//       with going to a different number of processors
//       (in collated format). Use file-based searching instead

namespace Foam
{

// If indeed collated format:
// Collect block-number in individual filenames
// (might differ on different processors)
static bool checkFileExistenceCollated
(
    const Foam::fileOperation& handler,
    const Foam::fileName& fName
)
{
    // using namespace Foam;

    bool found = false;
    {
        const label handlerComm = handler.comm();

        const label globalProci = UPstream::myProcNo(UPstream::worldComm);
        const label handlerProci = UPstream::myProcNo(handlerComm);
        const label nHandlerProcs = UPstream::nProcs(handlerComm);

        // Determine my local block number
        label myBlockNumber = -1;
        {
            fileOperation::procRangeType group;
            label proci = fileOperation::detectProcessorPath(fName, group);

            if (proci == -1 && group.empty())
            {
                // 'processorsXXX' format so contains all ranks
                // according to worldComm
                myBlockNumber = globalProci;
            }
            else
            {
                // 'processorsXXX_n-m' format so check for relative rank
                myBlockNumber = handlerProci;
            }
        }

        // Since we are streaming anyhow, could also pack as tuple:
        // Tuple2<fileName, label>


        // Collect file names on master of local communicator
        const fileNameList fNames
        (
            Pstream::listGatherValues
            (
                fName,
                handlerComm,
                UPstream::msgType()
            )
        );

        // Collect block numbers on master of local communicator
        const labelList myBlockNumbers
        (
            Pstream::listGatherValues
            (
                myBlockNumber,
                handlerComm,
                UPstream::msgType()
            )
        );

        // Determine for all whether the filename exists in the collated file.
        boolList allFound;

        if (UPstream::master(handlerComm))
        {
            allFound.resize(nHandlerProcs, false);

            // Store nBlocks and index of file that was used for nBlocks
            label nBlocks = -1;
            label blockRanki = -1;
            forAll(fNames, ranki)
            {
                if
                (
                    blockRanki == -1
                 || (fNames[ranki] != fNames[blockRanki])
                )
                {
                    blockRanki = ranki;
                    IFstream is(fNames[ranki]);
                    nBlocks = decomposedBlockData::getNumBlocks(is);
                }

                allFound[ranki] = (myBlockNumbers[ranki] < nBlocks);
            }
        }

        // Scatter using the handler communicator
        found = Pstream::listScatterValues
        (
            allFound,
            handlerComm,
            UPstream::msgType()
        );
    }

    return found;
}

} // End namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::bitSet Foam::haveProcessorFile
(
    const word& name,           // eg "faces"
    const fileName& instance,   // eg "constant"
    const fileName& local,      // eg, polyMesh
    const Time& runTime,
    const bool verbose
)
{
    const auto& handler = Foam::fileHandler();

    const fileName fName
    (
        handler.filePath(runTime.path()/instance/local/name)
    );

    bool found = handler.isFile(fName);

    // Assume non-collated (as fallback value).
    // If everyone claims to have the file, use master to verify if
    // collated is involved.

    bool isCollated = false;

    if (returnReduceAnd(found, UPstream::worldComm))
    {
        // Test for collated format.
        // Note: can test only world-master. Since even host-collated will have
        // same file format type for all processors
        if (UPstream::master(UPstream::worldComm))
        {
            const bool oldParRun = UPstream::parRun(false);

            if (IFstream is(fName); is.good())
            {
                IOobject io(name, instance, local, runTime);
                io.readHeader(is);

                isCollated = decomposedBlockData::isCollatedType(io);
            }

            UPstream::parRun(oldParRun);
        }
        Pstream::broadcast(isCollated, UPstream::worldComm);
    }

    // For collated, check that the corresponding blocks exist
    if (isCollated)
    {
        found = checkFileExistenceCollated(handler, fName);
    }


    // Globally consistent information about who has the file
    bitSet haveFileOnProc = bitSet::allGather(found, UPstream::worldComm);

    if (verbose)
    {
        Info<< "Per processor availability of \""
            << name << "\" file in " << instance/local << nl
            << "    " << flatOutput(haveFileOnProc) << nl << endl;
    }

    return haveFileOnProc;
}


Foam::boolList Foam::haveMeshFile
(
    const word& name,           // eg "faces"
    const fileName& instance,   // eg "constant"
    const fileName& local,      // eg, polyMesh
    const Time& runTime,
    const bool verbose
)
{
    const auto& handler = Foam::fileHandler();

    const fileName fName
    (
        handler.filePath(runTime.path()/instance/local/name)
    );

    bool found = handler.isFile(fName);

    // Assume non-collated (as fallback value).
    // If everyone claims to have the file, use master to verify if
    // collated is involved.

    bool isCollated = false;

    if (returnReduceAnd(found, UPstream::worldComm))
    {
        // Test for collated format.
        // Note: can test only world-master. Since even host-collated will have
        // same file format type for all processors
        if (UPstream::master(UPstream::worldComm))
        {
            const bool oldParRun = UPstream::parRun(false);

            if (IFstream is(fName); is.good())
            {
                IOobject io(name, instance, local, runTime);
                io.readHeader(is);

                isCollated = decomposedBlockData::isCollatedType(io);
            }

            UPstream::parRun(oldParRun);
        }
        Pstream::broadcast(isCollated, UPstream::worldComm);
    }

    // For collated, check that the corresponding blocks exist
    if (isCollated)
    {
        found = checkFileExistenceCollated(handler, fName);
    }


    // Globally consistent information about who has a mesh
    boolList haveFileOnProc
    (
        UPstream::allGatherValues<bool>(found, UPstream::worldComm)
    );

    if (verbose)
    {
        Info<< "Per processor availability of \""
            << name << "\" file in " << instance/local << nl
            << "    " << flatOutput(haveFileOnProc) << nl << endl;
    }

    return haveFileOnProc;
}


Foam::boolList Foam::haveMeshFile
(
    const Time& runTime,
    const fileName& meshPath,
    const word& meshFile,
    const bool verbose
)
{
    const auto& handler = Foam::fileHandler();

    const fileName fName
    (
        handler.filePath(runTime.path()/meshPath/meshFile)
    );

    bool found = handler.isFile(fName);

    // Assume non-collated (as fallback value).
    // If everyone claims to have the file, use master to verify if
    // collated is involved.

    bool isCollated = false;

    if (returnReduceAnd(found, UPstream::worldComm))
    {
        // Test for collated format.
        // Note: can test only world-master. Since even host-collated will have
        // same file format type for all processors
        if (UPstream::master(UPstream::worldComm))
        {
            const bool oldParRun = UPstream::parRun(false);

            if (IFstream is(fName); is.good())
            {
                IOobject io(meshFile, meshPath, runTime);
                io.readHeader(is);

                isCollated = decomposedBlockData::isCollatedType(io);
            }

            UPstream::parRun(oldParRun);
        }
        Pstream::broadcast(isCollated, UPstream::worldComm);
    }

    // For collated, check that the corresponding blocks exist
    if (isCollated)
    {
        found = checkFileExistenceCollated(handler, fName);
    }


    // Globally consistent information about who has a mesh
    boolList haveFileOnProc
    (
        UPstream::allGatherValues<bool>(found, UPstream::worldComm)
    );

    if (verbose)
    {
        Info<< "Per processor availability of \""
            << meshFile << "\" file in " << meshPath << nl
            << "    " << flatOutput(haveFileOnProc) << nl << endl;
    }

    return haveFileOnProc;
}


void Foam::removeProcAddressing(const faMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "edge", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::removeProcAddressing(const polyMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "cell", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::masterMeshInstance
(
    const IOobject& io,
    fileName& facesInstance,
    fileName& pointsInstance
)
{
    const fileName meshSubDir
    (
        polyMesh::meshDir(io.name())
    );

    if (UPstream::master())
    {
        const bool oldParRun = UPstream::parRun(false);
        const label oldNumProcs = fileHandler().nProcs();
        const int oldCache = fileOperation::cacheLevel(0);

        facesInstance = io.time().findInstance
        (
            meshSubDir,
            "faces",
            IOobjectOption::MUST_READ
        );
        pointsInstance = io.time().findInstance
        (
            meshSubDir,
            "points",
            IOobjectOption::MUST_READ
        );

        fileOperation::cacheLevel(oldCache);
        if (oldParRun)
        {
            const_cast<fileOperation&>(fileHandler()).nProcs(oldNumProcs);
        }
        UPstream::parRun(oldParRun);
    }

    // Broadcast information to all
    Pstream::broadcasts
    (
        UPstream::worldComm,
        facesInstance,
        pointsInstance
    );
}


// ************************************************************************* //
