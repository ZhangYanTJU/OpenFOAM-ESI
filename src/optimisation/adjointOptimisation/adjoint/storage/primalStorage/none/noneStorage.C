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

#include "noneStorage.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noneStorage, 0);
    defineRunTimeSelectionTable(noneStorage, dictionary);
    addToRunTimeSelectionTable(primalStorage, noneStorage, primalStorage);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noneStorage::noneStorage
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
:
    primalStorage(primalSolverObj, storageDict)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::noneStorage> Foam::noneStorage::New
(
    primalSolver& primalSolverObj,
    const dictionary storageDict
)
{
    return autoPtr<noneStorage>
    (
        new noneStorage(primalSolverObj, storageDict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::noneStorage::scratch()
{
    // Does nothing
}


void Foam::noneStorage::storeVariables()
{
    // Does nothing
}


void Foam::noneStorage::retrieveVariables()
{
    // Does nothing
}


// ************************************************************************* //
