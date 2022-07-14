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

#include "primalSolver.H"
#include "shortFullStorage.H"
#include "memInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shortFullStorage, 0);
    addToRunTimeSelectionTable
    (
        compressedFullStorage,
        shortFullStorage,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::shortFullStorage::gatherMetrics()
{
    if (storageParams_().timing())
    {
        for (scalar& storageMetric : storageMetrics_)
        {
            reduce(storageMetric, sumOp<scalar>());
        }
    }
}


void Foam::shortFullStorage::setMetricsFilesPtr()
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
            << setw(width) << "Initial Size (Mb)" << tab
            << setw(width) << "Initial Size Saved (Mb)" << endl;
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


inline void Foam::shortFullStorage::writeToFile()
{
    if (Pstream::master())
    {
        if (!metricsFilePtr_)
        {
            setMetricsFilesPtr();
        }
        const scalar initialSize = storageMetrics_[0];
        const scalar uncompressedRoughSize = storageMetrics_[1];
        unsigned int width = IOstream::defaultPrecision() + 6;
        metricsFilePtr_()
            << setprecision(IOstream::defaultPrecision());
        metricsFilePtr_() << setw(width) << initialSize*1.e-6 << tab
            << setw(width) << uncompressedRoughSize*1.e-6 << endl;
    }
}


void Foam::shortFullStorage::calcStorageMetrics()
{
    const scalarList& metrics = allPrimalVars_[iPtr_].storageMetrics();
    if (allPrimalVars_.size() == 1)
    {
        storageMetrics_ = metrics;
    }
    else
    {
        storageMetrics_ = storageMetrics_ + metrics;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortFullStorage::shortFullStorage
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    compressedFullStorage(primalSolverObj, storageDict),
    iPtr_(storagePtrs_[0])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shortFullStorage::scratch()
{
    totalSteps_ = 0;
    iPtr_ = -1;
    allPrimalVars_.clear();
    storageMetrics_.clear();
}


void Foam::shortFullStorage::storeVariables()
{
    if (allPrimalVars_.size() > iPtr_ + 1)
    {
        Info<< endl;
        WarningInFunction
            << tab
            << "allPrimalVars.size(), iPtr = "
            << allPrimalVars_.size() << " " << iPtr_
            << nl << tab
            << "Normally, allPrimalVars.size() = iPtr + 1 in each time step"
            << nl << tab
            << "Possible cause: Multiple compressions/decompressions "
            << "of the same time step"
            << nl << tab
            << "allPrimalVars.size() manually set to "
            << iPtr_ + 1 << nl << endl;
        allPrimalVars_.resize(iPtr_ + 1);
    }
    incompressibleVars& incoVars = refCast<incompressibleVars>(variablesSet_);
    autoPtr<compressedIncompressibleVars> varsPtr
    (
        compressedIncompressibleVars::New
        (
            incoVars,
            storageParams_()
        )
    );
    varsPtr->compress();
    varsPtr->calculateStorageMetrics();
    allPrimalVars_.append(varsPtr);
    ++iPtr_;
    calcStorageMetrics();
    // Accumulate the primal step to counter
    totalSteps_++;
}


void Foam::shortFullStorage::retrieveVariables()
{
    if (allPrimalVars_.empty())
    {
        FatalErrorInFunction
            << "allPrimalVars_ list is empty." << endl
            << exit(FatalError);
    }
    allPrimalVars_[iPtr_].decompress();
    variablesSet_.validateTurbulence();
    --iPtr_;
    // Accumulate the adjoint step to counter
    totalSteps_++;
}


void Foam::shortFullStorage::storageMetrics()
{
    gatherMetrics();
    writeToFile();

    const scalar initialSize = storageMetrics_[0];
    const scalar uncompressedRoughSize = storageMetrics_[1];

    Info<< nl << "- - - - - - - - - - - - - - - - - - - - -" << endl;
    Info<< "Primal Storage Statistics" << endl;
    Info<< "- - - - - - - - - - - - - - - - - - - - -" << nl << endl;

    Info<< setprecision(6);
    Info<< "Initial Size = " << initialSize*1.e-6 << " Mb" << endl;
    Info<< "Stored Size  = " << uncompressedRoughSize*1.e-6  << " Mb" << endl;
    Info<< setprecision(IOstream::defaultPrecision()) << endl;
}


// ************************************************************************* //
