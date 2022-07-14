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

#include "Ostream.H"
#include "binomialCheckPointing.H"
#include "memInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "messageStream.H"
#include "primalSolver.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binomialCheckPointing, 0);
    defineRunTimeSelectionTable(binomialCheckPointing, dictionary);
    addToRunTimeSelectionTable
    (
        primalStorage,
        binomialCheckPointing,
        primalStorage
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::binomialCheckPointing::hostStartProcess
(
    wordList& hostNames
) const
{
    // Gather host names for all processors
    wordList machineName(Pstream::nProcs());
    machineName[Pstream::myProcNo()] = hostName();
    Pstream::gatherList(machineName);
    Pstream::scatterList(machineName);
    hostNames = machineName;

    // Start processor ID for each host; to be resized
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
            iHost++;
        }
    }
    // Allocate an additional entry for easy looping over start
    // of this host and the next one
    hostStartProcess.setSize(iHost + 1);
    hostStartProcess[iHost] = nProcs;
    return hostStartProcess;
}


Foam::scalar Foam::binomialCheckPointing::checkPointSize()
{
    Time& time = const_cast<Time&>(mesh_.time());
    label startTimeIndex = mesh_.time().timeIndex();
    dimensionedScalar startTime = mesh_.time();

    Info<< "CheckPoint size computation:: Time = "
        << time.timeName() << nl << endl;
    primalSolver_.solveIter();
    checkPoint temp(storageParams_());
    incompressibleVars& incoVars = refCast<incompressibleVars>(variablesSet_);
    temp.setPlaceHolder(false);
    temp.store(incoVars);

    time.setTime(startTime, startTimeIndex);

    return temp.storageMetrics()[1];
}


unsigned long long Foam::binomialCheckPointing::maxNumberOfCheckPoints()
{
    // Find checkpoint size per processor (in MB)
    scalarList size(Pstream::nProcs());
    size[Pstream::myProcNo()] = checkPointSize()*1.e-6;
    Pstream::gatherList(size);
    Pstream::scatterList(size);
    if (!memInfo().valid())
    {
        FatalErrorInFunction
            << "memInfo is invalid" << endl
            << exit(FatalError);
    }
    // Available memory in MB. Is stored on a per-processor base but actually
    // refers to the host. Use with caution
    scalarList free(Pstream::nProcs());
    free[Pstream::myProcNo()] = memInfo().update().free()*1.e-3;
    Pstream::gatherList(free);
    Pstream::scatterList(free);

    // Host name, per process
    wordList hostNames(0);
    // Start processor ID per host
    labelList hostStartProcessID(hostStartProcess(hostNames));
    label nHosts(hostStartProcessID.size() - 1);

    // Determine the number of checkpoints that can be allocated as the
    // minimum from all hosts
    labelList hostCheckPointSize(nHosts, 0);
    label nCheckPoints(pTraits<label>::max);
    label limitingHostID(-1);
    for (label iHost = 0; iHost < nHosts; ++iHost)
    {
        const label start(hostStartProcessID[iHost]);
        const label end(hostStartProcessID[iHost + 1]);
        label& hostSize = hostCheckPointSize[iHost];
        for (label procI = start; procI < end; ++procI)
        {
            hostSize += size[procI];
        }
        if (hostSize == 0)
        {
            FatalErrorInFunction
                << "hostSize = 0!!" << endl
                << exit(FatalError);
        }
        scalar val = free[start]/hostSize;
        label nHostCheckPoints(pTraits<label>::max);
        if (val <= 0)
        {
            FatalErrorInFunction
                << "free[start], hostSize = " << free[start]
                << " " << hostSize << endl
                << exit(FatalError);
        }
        else if (val < pTraits<label>::max)
        {
            nHostCheckPoints = val;
        }
        if (nHostCheckPoints < nCheckPoints)
        {
            nCheckPoints = nHostCheckPoints;
            limitingHostID = iHost;
        }
    }

    // Print global check point size
    Info<< nl << "Global size of a checkpoint is "
        << sum(size) << " MBs" << nl
        << "Total free memory in all hosts: "
        << sum(free)*1.e-3 << " GBs" << endl;

    // Print max number of check points that fit in the most limiting host
    Info<< "Approximately " << nCheckPoints
        << " checkpoints can be stored, "
        << "without accounting for cashed memory " << endl;

    // Print information on a per host basis
    if (debug > 1)
    {
        Info<< nl
            << setw(23) << "Host" << " | "
            << setw(23) << "Checkpoint size (kBs)" << " | "
            << setw(23) << "Available memory (kBs)" << " | "
            << setw(23) << "Max checkpoints (local)"
            << nl;
        for (label iHost = 0; iHost < nHosts; ++iHost)
        {
            const label start(hostStartProcessID[iHost]);
            Info<< setw(23) << hostNames[start] << "   "
                << setw(23) << hostCheckPointSize[iHost] << "   "
                << setw(23) << free[start] << "   "
                << setw(23) << free[start]/hostCheckPointSize[iHost]
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

    return nCheckPoints;
}


void Foam::binomialCheckPointing::initialize()
{
    nCheckPoints_ = storageDict_.subDict("storageParameters").
        get<label>("nCheckPoints") + 1;
    if (nCheckPoints_ < 1)
    {
        FatalErrorInFunction
            << "nCheckPoints = " << nCheckPoints_ << ". Please provide a "
            << "value >=1" << endl
            << exit(FatalError);
    }
    /*if (storageDict_.subDict("storageParameters").found("nCheckPoints"))
    {
        nCheckPoints_ = storageDict_.subDict("storageParameters").
            get<label>("nCheckPoints");
    }
    else
    {
        nCheckPoints_ = maxNumberOfCheckPoints();
    }*/
    checkPoints_.setSize(nCheckPoints_);
    indirectAddressing_ = identity(nCheckPoints_);
    Info<< nl << "Total number of checkpoints = "
        << nCheckPoints_ << nl << endl;
    for (label fI = 0; fI < nCheckPoints_; ++fI)
    {
        checkPoints_.set(fI, new checkPoint(storageParams_()));
    }
}


void Foam::binomialCheckPointing::gatherMetrics()
{
    if (storageParams_().timing())
    {
        for (scalar& smI : storageMetrics_)
        {
            reduce(smI, sumOp<scalar>());
        }
    }
}


void Foam::binomialCheckPointing::setMetricsFilesPtr()
{
    fileName filePath =
        storageFolder_/"AllOptCycles" + primalSolver_.solverName();
    if (!isFile(filePath))
    {
        DebugInfo
            << "Creating file: " << filePath << endl;
        metricsFilePtr_.reset(new OFstream(filePath));
        if (!metricsFilePtr_()())
        {
            FatalErrorInFunction
                << "Stream has failed." << endl
                << exit(FatalError);
        }
        if (!isFile(filePath))
        {
            FatalErrorInFunction
                << "Error while creating file: " << filePath << endl
                << exit(FatalError);
        }
        unsigned int width = IOstream::defaultPrecision() + 6;
        metricsFilePtr_()
            << setw(width) << "nCheckPoints" << tab
            << setw(width) << "Overall CR (CR_0)" << tab
            << setw(width) << "ZFP CR (CR_ZFP)" << tab
            << setw(width) << "Initial Size (Mb)" << tab
            << setw(width) << "Initial Size Saved (Mb)" << tab
            << setw(width) << "Compressed Size (Mb)"<< endl;
    }
    else
    {
        metricsFilePtr_.reset
        (
            new OFstream
                (
                    filePath,
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::UNCOMPRESSED,
                    true
                )
        );
    }
}


void Foam::binomialCheckPointing::setCheckPointsFilesPtr()
{
    fileName filePath =
        storageFolder_/"checkPoints" + primalSolver_.solverName();
    if (!isFile(filePath))
    {
        DebugInfo
            << "Creating file: " << filePath << endl;
        checkPointsFilePtr_.reset(new OFstream(filePath));
        if (!checkPointsFilePtr_()())
        {
            FatalErrorInFunction
                << "Stream has failed." << endl
                << exit(FatalError);
        }
        if (!isFile(filePath))
        {
            FatalErrorInFunction
                << "Error while creating file: " << filePath << endl
                << exit(FatalError);
        }
        unsigned int width = IOstream::defaultPrecision() + 6;
        checkPointsFilePtr_()
            << setw(width) << "S/N" << tab
            << setw(width) << "TimeIndex" << tab
            << setw(width) << "Level" << tab
            << setw(width) << "Time" << tab
            << setw(width) << "Placeholder" << endl;
    }
    else
    {
        checkPointsFilePtr_.reset
        (
            new OFstream
                (
                    filePath,
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::UNCOMPRESSED,
                    true
                )
        );
    }
}


inline void Foam::binomialCheckPointing::writeToFile()
{
    if (Pstream::master())
    {
        if (!metricsFilePtr_)
        {
            setMetricsFilesPtr();
        }
        const scalar& initialSize = storageMetrics_[0];
        const scalar& uncompressedRoughSize = storageMetrics_[1];
        const scalar& uncompressedSize = storageMetrics_[2];
        const scalar& compressedSize = storageMetrics_[3];

        unsigned int width = IOstream::defaultPrecision() + 6;
        metricsFilePtr_()
            << setprecision(IOstream::defaultPrecision());
        metricsFilePtr_()
            << setw(width) << activeCheckPoints_ << tab
            << setw(width) << uncompressedRoughSize/compressedSize << tab
            << setw(width) << uncompressedSize/compressedSize << tab
            << setw(width) << initialSize*1.e-6 << tab
            << setw(width) << uncompressedRoughSize*1.e-6 << tab
            << setw(width) << compressedSize*1e-6 << endl;
        if (!checkPointsFilePtr_)
        {
            setCheckPointsFilesPtr();
        }
        checkPointsFilePtr_()
            << setprecision(IOstream::defaultPrecision());
        for (label fI = 0; fI < nCheckPoints_; ++fI)
        {
            const label pos = indirectAddressing_[fI];
            if (checkPoints_[pos].active())
            {
                checkPointsFilePtr_()
                    << setw(width) << fI << tab
                    << setw(width) << checkPoints_[pos].timeIndex() << tab
                    << setw(width) << checkPoints_[pos].level() << tab
                    << setw(width) << name(checkPoints_[pos].timeValue())
                    << setw(width) << checkPoints_[pos].placeHolder()
                    << endl;
            }
        }
        checkPointsFilePtr_() << endl;
    }
}


void Foam::binomialCheckPointing::calcStorageMetrics()
{
    storageMetrics_ = checkPoints_[0].storageMetrics();
    for (label fI = 1; fI < nCheckPoints_; ++fI)
    {
        if (checkPoints_[fI].active() && !checkPoints_[fI].placeHolder())
        {
            storageMetrics_ =
                storageMetrics_
              + checkPoints_[fI].storageMetrics();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binomialCheckPointing::binomialCheckPointing
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    primalStorage(primalSolverObj, storageDict),
    totalSteps_(0),
    recomputationSteps_(0),
    storageParams_
    (
        storageParameters::New
            (mesh_, storageDict_, variablesSet_.allocatedFieldNames())
    ),
    adjustTimeStep_(storageParams_->adjustTimeStep()),
    checkPoints_(),
    startingTime_(mesh_.time().startTime()),
    activeCheckPoints_(0),
    iPtr_(storagePtrs_[0])
{
    makeFolder();
    initialize();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::binomialCheckPointing> Foam::binomialCheckPointing::New
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
{
    return autoPtr<binomialCheckPointing>
    (
        new binomialCheckPointing(primalSolverObj, storageDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::binomialCheckPointing::scratch()
{
    indirectAddressing_ = identity(nCheckPoints_);
    for (checkPoint& cp : checkPoints_)
    {
        cp.empty();
    }
    totalSteps_ = 0;
    recomputationSteps_ = 0;
    activeCheckPoints_ = 0;
    storageMetrics_.clear();
    iPtr_ = -1;
}


void Foam::binomialCheckPointing::storeVariables()
{
    label jPtr(0);
    // Add current to placeholder

    // Find the dispensable checkpoint with the largest time index
    // idis holds the index of the checkpoint that will be removed.
    label maxlevel = 0;
    // Despensable checkpoint index
    label idis = 0;
    for (label i = nCheckPoints_ - 1; i >= 1; --i)
    {
        if (maxlevel > checkPoints_[indirectAddressing_[i]].level())
        {
            // Keep the index only the first time it enters, since it
            // will have the greater index
            idis = i;
            break;
        }
        else
        {
            maxlevel = checkPoints_[indirectAddressing_[i]].level();
        }
    }
    if (activeCheckPoints_ < nCheckPoints_)
    {
        // Scenario 1: not all checkpoints allocated yet.
        DebugInfo
            << "Scenario 1" << endl;
        const label pos = indirectAddressing_[activeCheckPoints_];
        checkPoints_[pos].setPlaceHolder(false);
        jPtr = activeCheckPoints_;
        activeCheckPoints_++;
    }
    else if (idis > 0)
    {
        // Scenario 2:
        // at least one checkpoint is despensable. Create a new checkPoint
        // at the end of the indirect addressing list of level 0.
        DebugInfo
            << "Scenario 2" << endl;
        const label pos = indirectAddressing_[idis];
        for (label i = idis; i < nCheckPoints_ - 1; ++i)
        {
            indirectAddressing_[i] = indirectAddressing_[i + 1];
        }
        indirectAddressing_[nCheckPoints_ - 1] = pos;
        checkPoints_[pos].setPlaceHolder(false);
        jPtr = nCheckPoints_ - 1;
    }
    else
    {
        // Scenario 3:
        // no checkpoint is despensable. Overwrite last checkpoint
        // and increase its level
        DebugInfo
            << "Scenario 3" << endl;
        const label pos = indirectAddressing_[nCheckPoints_ - 1];
        checkPoints_[pos].setPlaceHolder(true);
        jPtr = nCheckPoints_ - 1;
    }

    // If the previous time-step is in the current set of check-points and
    // checkpoint is not a placeholder checkpoint store the solution. When
    // retrieving the solution at scenario 2, the checkpoint is not set as a
    // placeholder as Wang proposes, so there is no need to store again the
    // solution in that case.
    if
    (
        (
            time_.timeIndex() - 1
         == checkPoints_[indirectAddressing_[jPtr - 1]].timeIndex()
        )
     && checkPoints_[indirectAddressing_[jPtr - 1]].placeHolder()
    )
    {
        DebugInfo
            << "Storing time-step with time, timeIndex = "
            << name(checkPoints_[indirectAddressing_[jPtr - 1]].timeValue())
            << " "
            << checkPoints_[indirectAddressing_[jPtr - 1]].timeIndex()
            << endl;
        incompressibleVars& incoVars =
            refCast<incompressibleVars>(variablesSet_);
        checkPoints_[indirectAddressing_[jPtr - 1]].store(incoVars);
    }
    calcStorageMetrics();
    ++iPtr_;
    if (debug)
    {
        Info<< "checkPoint time values -->Start<--" << endl;
        for (label fI = 0; fI < nCheckPoints_; ++fI)
        {
            const label pos = indirectAddressing_[fI];
            Info<< "Checkpoint-ID, indirectAddressing, timeIndex, "
                << "level, time, placeHolder = "
                << fI << " "
                << pos << " "
                << checkPoints_[pos].timeIndex() << " "
                << checkPoints_[pos].level() << " "
                << name(checkPoints_[pos].timeValue()) << " "
                << checkPoints_[pos].placeHolder() << endl;
        }
        Info<< "checkPoint time values --> End <--" << endl;
    }
    // Accumulate the primal step to counter
    totalSteps_++;
}


void Foam::binomialCheckPointing::storeInitialVariables()
{
    // Store initial time-step as a placeholder checkPoint of level oo
    checkPoints_[0].setPlaceHolder(false);
    checkPoints_[0].level() = pTraits<label>::max;
    activeCheckPoints_++;
    ++iPtr_;
}


void Foam::binomialCheckPointing::retrieveVariables()
{
    DebugInfo
        << "Retrieving checkPoint for time " << time_.timeIndex() << endl;
    dimensionedScalar adjEndingTime = time_.endTime();
    label timeIndex = time_.timeIndex();

    // Remove the checkpoint if it is at a greater time-step than the current
    // one, or if it is at the same time-step as the current one, but it is a
    // placeHolder checkPoint.
    const label lastPos = indirectAddressing_[activeCheckPoints_ - 1];
    if
    (
        (checkPoints_[lastPos].timeIndex() > timeIndex)
     || (
            checkPoints_[lastPos].timeIndex() == timeIndex &&
            checkPoints_[lastPos].placeHolder()
        )
    )
    {
        const label pos = indirectAddressing_[activeCheckPoints_ - 1];
        DebugInfo
            << "Removing checkPoint: ID, timeIndex, time = " << pos << " "
            << checkPoints_[pos].timeIndex() << " "
            << name(checkPoints_[pos].timeValue()) << endl;
        checkPoints_[pos].empty();
        activeCheckPoints_--;
    }
    Time& time = const_cast<Time&>(time_);
    if
    (
        checkPoints_[indirectAddressing_[activeCheckPoints_ - 1]].timeIndex()
     == timeIndex
    )
    {
        // Scenario 1:
        // there is a checkpoint at the current time-step. Retrieve the
        // solution and make it a placeholder checkpoint
        DebugInfo
            << "Scenario 1" << endl;
        const label pos = indirectAddressing_[activeCheckPoints_ - 1];
        checkPoints_[pos].retrieve();
        variablesSet_.validateTurbulence();
        if (adjustTimeStep_)
        {
            time.setDeltaT(checkPoints_[pos].deltaT());
        }
        checkPoints_[pos].resetToPlaceHolder();
        activeCheckPoints_--;
    }
    else
    {
        // Scenario 2:
        // there is no checkpoint at the current time-step. Retrieve the
        // solution at the last checkpoint
        DebugInfo
            << "Scenario 2" << endl;
        const label pos = indirectAddressing_[activeCheckPoints_ - 1];

        scalar nowTime(time.value());
        if (adjustTimeStep_)
        {
            time.setDeltaT(checkPoints_[pos].deltaT());
        }
        time.setTime
            (checkPoints_[pos].timeValue(), checkPoints_[pos].timeIndex());
        time.setEndTime(nowTime);

        checkPoints_[pos].retrieve();
        variablesSet_.validateTurbulence();

        while (primalSolver_.loop())
        {
            Info<< "Primal-CheckPointing::";
            primalSolver_.solveIter();
            recomputationSteps_++;
        }
    }
    time.setEndTime(adjEndingTime);
    --iPtr_;
    // Accumulate the adjoint step to counter
    totalSteps_++;
}


void Foam::binomialCheckPointing::postAdjointLoop()
{
    // Output re-computation info
    scalar recomputationLoad =
        scalar((2.*recomputationSteps_))/
        scalar((totalSteps_-recomputationSteps_));
    Info<<"\n/* * * * * * * * * * * * * * * * * * * * * * * */" << nl
        <<"Total number of steps "<<totalSteps_ << nl
        <<"Recomputations        "<<recomputationSteps_ << nl
        <<"Recomputation load    "<<recomputationLoad << nl
        <<"/* * * * * * * * * * * * * * * * * * * * * * * */" << nl << endl;

    // Rewind storage pointers with a Warning
    WarningInFunction
        << "Rewinding storage pointers, "
        << "even though check-points have been re-arranged." << nl
        << "Use binomialCheckPointing with caution in the presence of multiple"
        << " adjoint solvers"
        << endl;
    primalStorage::postAdjointLoop();
}


void Foam::binomialCheckPointing::storageMetrics()
{
    gatherMetrics();
    writeToFile();

    const scalar initialSize = storageMetrics_[0];
    const scalar uncompressedRoughSize = storageMetrics_[1];
    const scalar uncompressedSize = storageMetrics_[2];
    const scalar compressedSize = storageMetrics_[3];

    Info<< nl << "- - - - - - - - - - - - - - - - - - - - -" << endl;
    Info<< "Primal Storage Statistics" << endl;
    Info<< "- - - - - - - - - - - - - - - - - - - - -" << nl << endl;

    unsigned int width = 6;//IOstream::defaultPrecision();
    Info<< setprecision(width);
    Info<< "nCheckPoints                  = "
        << activeCheckPoints_ << nl
        << "Compression Ratio CR_0       = "
        << uncompressedRoughSize/compressedSize << nl
        << "ZFP Compression ratio CR_ZFP = "
        << uncompressedSize/compressedSize << nl
        << "Initial Size                 = "
        << initialSize*1.e-6 << " Mb" << nl
        << "Initial Size Saved           = "
        << uncompressedRoughSize*1.e-6  << " Mb" << nl
        << "Compressed Size              = "
        << compressedSize*1.e-6 << " Mb" << endl;
    Info<< setprecision(IOstream::defaultPrecision()) << endl;
    /*
    Info << "checkPoint time values -->Start<--" << endl;
    for (label fI=0; fI<nCheckPoints_; fI++)
    {
        const label& pos = indirectAddressing_[fI];
        Info << "Checkpoint-ID, timeIndex, level, time, placeholder = " << fI << " "
             << checkPoints_[pos].timeIndex() << " " << checkPoints_[pos].level()
             << " " << ::Foam::name(checkPoints_[pos].timeValue()) << " "
             << checkPoints_[pos].placeHolder() << endl;
    }
    Info << "checkPoint time values --> End <--" << endl;
    */
}


// ************************************************************************* //
