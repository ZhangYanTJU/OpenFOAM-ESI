/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(primalSolver, 0);
    defineRunTimeSelectionTable(primalSolver, primalSolver);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::primalSolver::primalSolver
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    Switch useCustomReadTime
)
:
    solver(mesh, managerType, dict),
    primalStorage_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::primalSolver> Foam::primalSolver::New
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    Switch useCustomReadTime
)
{
    const word solverType(dict.get<word>("type"));

    auto* ctorPtr = primalSolverConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "primalSolver",
            solverType,
            *primalSolverConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return
        autoPtr<primalSolver>
        (
            ctorPtr(mesh, managerType, dict, useCustomReadTime)
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primalSolver::readDict(const dictionary& dict)
{
    if (solver::readDict(dict))
    {
        return true;
    }

    return false;
}


Foam::autoPtr<Foam::primalStorage>& Foam::primalSolver::getPrimalStorage()
{
    return primalStorage_;
}


void Foam::primalSolver::correctBoundaryConditions()
{
    // Do nothing
}


// ************************************************************************* //
