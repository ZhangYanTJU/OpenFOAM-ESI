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

#include "error.H"
#include "primalStorage.H"
#include "primalSolver.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(primalStorage, 0);
    defineRunTimeSelectionTable(primalStorage, primalStorage);

    label primalStorage::counter = 0;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::primalStorage::makeFolder()
{
    if (Pstream::master())
    {
        storageFolder_ = time_.globalPath()/"optimisation"/"storage";
        //- do not remove directory to allow for continuation
        //if (isDir(storageFolder_)) rmDir(storageFolder_);
        mkDir(storageFolder_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primalStorage::primalStorage
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    time_(primalSolverObj.mesh().time()),
    mesh_(primalSolverObj.mesh()),
    primalSolver_(primalSolverObj),
    storageDict_(storageDict),
    variablesSet_(primalSolver_.getVariablesSet()),
    storagePtrs_(1, -1),
    storagePtrsCopy_(nullptr)
{
    if (Pstream::master())
    {
        counter++;
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::primalStorage> Foam::primalStorage::New
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
{
    const word type(storageDict.get<word>("type"));
    Info<< "Storage type for the primal fields: " << type << endl;

    auto* ctorPtr = primalStorageConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            storageDict,
            "primalStorage",
            type,
            *primalStorageConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<primalStorage>(ctorPtr(primalSolverObj, storageDict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::primalStorage::storagePtrs() const
{
    return storagePtrs_;
}


void Foam::primalStorage::copyStoragePtrs()
{
    storagePtrsCopy_.reset(new labelList(storagePtrs_));
}


void Foam::primalStorage::rewindStoragePtrs()
{
    if (storagePtrsCopy_)
    {
        storagePtrs_ = storagePtrsCopy_();
    }
    else
    {
        FatalErrorInFunction
            << "Attempted to retrieve unset storePtrs" << endl;
    }
}


void Foam::primalStorage::storeInitialVariables()
{
    // Does nothing in base
}


void Foam::primalStorage::preLoopRetrieveVariables()
{
    // Does nothing in base
}


void Foam::primalStorage::preAdjointLoop()
{
    // Copy pointers to storage instances to a backup variable
    copyStoragePtrs();
    // Retrive certain variables before the adjoint loop
    // (e.g. from files, for restarted cases)
    preLoopRetrieveVariables();
}


void Foam::primalStorage::postAdjointLoop()
{
    rewindStoragePtrs();
}


void Foam::primalStorage::storageMetrics()
{
    // Does nothing in base
}


// ************************************************************************* //
