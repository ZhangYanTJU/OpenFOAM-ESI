/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "steadyOptimisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(steadyOptimisation, 0);
    addToRunTimeSelectionTable
    (
        optimisationManager,
        steadyOptimisation,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::steadyOptimisation::updateOptTypeSource()
{
    forAll(primalSolvers_, pI)
    {
        primalSolvers_[pI].updateOptTypeSource(optType_->sourcePtr());
    }

    forAll(adjointSolverManagers_, asmI)
    {
        PtrList<adjointSolver>& adjointSolvers =
            adjointSolverManagers_[asmI].adjointSolvers();

        forAll(adjointSolvers, aI)
        {
            adjointSolvers[aI].updateOptTypeSource(optType_->sourcePtr());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::steadyOptimisation::steadyOptimisation(fvMesh& mesh)
:
    optimisationManager(mesh)
{
    initialize();
    optType_.reset
    (
        incompressible::optimisationType::New
        (
            mesh,
            subDict("optimisation"),
            adjointSolverManagers_
        ).ptr()
    );

    // Update source ptrs in all solvers to look at the source held in optType
    // Possible problem if mesh is adapted
    updateOptTypeSource();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::optimisationManager& Foam::steadyOptimisation::operator++()
{
    time_++;
    if (!end())
    {
        Info<< "\n* * * * * * * * * * * * * * * * *" << endl;
        Info<< "Optimisation cycle " << time_.value() << endl;
        Info<< "* * * * * * * * * * * * * * * * *\n" << endl;
    }
    return *this;
}


Foam::optimisationManager& Foam::steadyOptimisation::operator++(int)
{
    return operator++();
}


bool Foam::steadyOptimisation::checkEndOfLoopAndUpdate()
{
    if (update())
    {
        optType_->update();
    }
    return end();
}


bool Foam::steadyOptimisation::end()
{
    return time_.end();
}


bool Foam::steadyOptimisation::update()
{
    return (time_.timeIndex() != 1 && !end());
}


// ************************************************************************* //
