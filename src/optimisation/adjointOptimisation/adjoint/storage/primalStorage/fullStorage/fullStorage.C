/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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

#include "fullStorage.H"
#include "memInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fullStorage, 0);
    defineRunTimeSelectionTable(fullStorage, dictionary);
    addToRunTimeSelectionTable(primalStorage, fullStorage, primalStorage);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fullStorage::hostStartProcess(wordList& hostNames) const
{
    // Gather host names for all processors
    wordList machineName(Pstream::nProcs());
    machineName[Pstream::myProcNo()] = hostName();
    Pstream::gatherList(machineName);
    Pstream::scatterList(machineName);
    hostNames = machineName;

    // Start process ID for each host; to be resized
    const label nProcs(Pstream::nProcs());
    labelList hostStartProcess(nProcs + 1, 0);

    label iHost(0);
    word prev(word::null);
    for (label procI = 0; procI < nProcs; ++procI)
    {
        if (machineName[procI] != prev)
        {
            hostStartProcess[iHost] = procI;
            prev = machineName[procI];
            ++iHost;
        }
    }
    // Allocate an additional entry for easy looping over start
    // of this host and the next one
    hostStartProcess.setSize(iHost + 1);
    hostStartProcess[iHost] = nProcs;
    return hostStartProcess;
}


Foam::label Foam::fullStorage::maxNumberOfInstances() const
{
    // Find instance size per processor
    labelList sizeBef(Pstream::nProcs());
    sizeBef[Pstream::myProcNo()] = memInfo().size();
    Pstream::gatherList(sizeBef);
    Pstream::scatterList(sizeBef);

    autoPtr<variablesSet> temp(variablesSet_.clone());

    labelList sizeAft(Pstream::nProcs());
    sizeAft[Pstream::myProcNo()] = memInfo().size();
    Pstream::gatherList(sizeAft);
    Pstream::scatterList(sizeAft);

    // Available memory. Is stored on a per-processor base but actually
    // refers to the host. Use with caution
    labelList free(Pstream::nProcs());
    free[Pstream::myProcNo()] = memInfo().free();
    Pstream::gatherList(free);
    Pstream::scatterList(free);

    // Host name, per processor basis
    wordList hostNames(0);
    // Start process ID per host
    labelList hostStartProcessID(hostStartProcess(hostNames));
    label nHosts(hostStartProcessID.size() - 1);

    // Determine the number of instances that can be allocated as the
    // minimum from all hosts
    labelList hostInstanceSize(nHosts, 0);
    label nInstances(pTraits<label>::max);
    label limitingHostID(-1);
    for (label iHost = 0; iHost < nHosts; ++iHost)
    {
        const label start(hostStartProcessID[iHost]);
        const label end(hostStartProcessID[iHost + 1]);
        label& hostSize = hostInstanceSize[iHost];
        for (label procI = start; procI < end; ++procI)
        {
            hostSize += sizeAft[procI] - sizeBef[procI];
        }
        label nHostCheckPoints(free[start]/hostSize);
        if (nHostCheckPoints < nInstances)
        {
            nInstances = nHostCheckPoints;
            limitingHostID = iHost;
        }
    }

    // Print global instance size
    const label checkPointSize(sum(sizeAft) - sum(sizeBef));
    Info<< nl << "Global size of an instance is "
        << checkPointSize << " KBs" << endl;

    // Print max number of check points that fit in the most limiting host
    Info<< "Approximately " << nInstances
        << " instances can be stored, "
        << "without accounting for cached memory " << endl;

    // Print information on a per host basis
    if (debug > 1)
    {
        Info<< nl
            << setw(23) << "Host" << " | "
            << setw(23) << "Instance size (kBs)" << " | "
            << setw(23) << "Available memory (kBs)" << " | "
            << setw(23) << "Max checkpoints (local)"
            << nl;
        for (label iHost = 0; iHost < nHosts; ++iHost)
        {
            const label start(hostStartProcessID[iHost]);
            Info<< setw(23) << hostNames[start] << "   "
                << setw(23) << hostInstanceSize[iHost] << "   "
                << setw(23) << free[start] << "   "
                << setw(23) << free[start]/hostInstanceSize[iHost]
                << endl;
        }
        Info<< endl;
    }

    string limitingHost;
    if (Pstream::myProcNo() == hostStartProcessID[limitingHostID])
    {
        limitingHost = hostName();
    }
    reduce(limitingHost, sumOp<word>());
    Info<< "Limiting factor is host " << limitingHost << nl << endl;

    return nInstances;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fullStorage::fullStorage
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    primalStorage(primalSolverObj, storageDict),
    allPrimalVars_(0),
    iPtr_(storagePtrs_[0])
{
    /*if (debug)
    {
        maxNumberOfCheckPoints();
    }*/
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fullStorage> Foam::fullStorage::New
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
{
    return autoPtr<fullStorage>
    (
        new fullStorage(primalSolverObj, storageDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fullStorage::scratch()
{
    allPrimalVars_.clear();
    iPtr_ = -1;
}


void Foam::fullStorage::storeVariables()
{
    allPrimalVars_.append(variablesSet_.clone());
    ++iPtr_;
}


void Foam::fullStorage::retrieveVariables()
{
    if (allPrimalVars_.empty())
    {
        FatalErrorInFunction
            << "allPrimalVars list is empty." << endl
            << exit(FatalError);
    }
    variablesSet_.transfer(allPrimalVars_[iPtr_]);
    //allPrimalVars_.release(iPtr_);
    --iPtr_;
}


// ************************************************************************* //
