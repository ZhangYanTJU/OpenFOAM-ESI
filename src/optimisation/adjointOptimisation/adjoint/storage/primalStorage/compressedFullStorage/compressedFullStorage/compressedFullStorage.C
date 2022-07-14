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

#include "compressedFullStorage.H"
#include "compressedIncompressibleVars.H"
#include "memInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressedFullStorage, 0);
    defineRunTimeSelectionTable(compressedFullStorage, dictionary);
    addToRunTimeSelectionTable
    (
        primalStorage,
        compressedFullStorage,
        primalStorage
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressedFullStorage::compressedFullStorage
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    primalStorage(primalSolverObj, storageDict),
    allPrimalVars_(0),
    storageMetrics_(),
    totalSteps_(0),
    recomputationSteps_(0),
    metricsFilePtr_()
{
    makeFolder();
    storageParams_.reset
    (
        storageParameters::New
        (
            mesh_,
            storageDict_,
            variablesSet_.allocatedFieldNames()
        )
    );
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::compressedFullStorage> Foam::compressedFullStorage::New
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
{
    const dictionary& dict = storageDict.subDict("storageParameters");
    const word type(dict.get<word>("algorithm") + "FullStorage");
    Info<< "CompressedFullStorage type" << tab << type << endl;

    auto* ctorPtr = dictionaryConstructorTable(type);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            storageDict,
            "compressedFullStorage",
            type,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<compressedFullStorage>
        (ctorPtr(primalSolverObj, storageDict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressedFullStorage::scratch()
{
    NotImplemented
}


void Foam::compressedFullStorage::storeVariables()
{
    NotImplemented
}


void Foam::compressedFullStorage::retrieveVariables()
{
    NotImplemented
}


void Foam::compressedFullStorage::postAdjointLoop()
{
    scalar recomputationLoad =
        scalar((2.*recomputationSteps_))/
        scalar((totalSteps_-recomputationSteps_));
    Info<<"\n/* * * * * * * * * * * * * * * * * * * * * * * */\n"
        <<"Total number of steps "<<totalSteps_ << nl
        <<"Recomputations        "<<recomputationSteps_ << nl
        <<"Recomputation load    "<<recomputationLoad << nl
        <<"/* * * * * * * * * * * * * * * * * * * * * * * */\n" << endl;

    // Rewind storage ptrs
    primalStorage::postAdjointLoop();
}


void Foam::compressedFullStorage::storageMetrics()
{
    NotImplemented
}


// ************************************************************************* //
