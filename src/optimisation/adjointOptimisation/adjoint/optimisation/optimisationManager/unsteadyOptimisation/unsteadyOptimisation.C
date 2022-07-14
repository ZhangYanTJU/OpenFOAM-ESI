/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "unsteadyOptimisation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(unsteadyOptimisation, 0);
    addToRunTimeSelectionTable
    (
        optimisationManager,
        unsteadyOptimisation,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::unsteadyOptimisation::updateOptTypeSource()
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::unsteadyOptimisation::resetTime()
{
    timeManip_.moveToPrimalStartTime();
}


void Foam::unsteadyOptimisation::moveDesignVariables()
{
    if (shouldUpdateDesignVariables_)
    {
        // Updated grid should be written at the primal end time, to coinside
        // with the primal and adjoint solutions to be stored then.  The
        // (unknown) time index is not important at this point, so it remains
        // unchanged.
        timeManip_.storeTime();
        timeManip_.moveToAdjointStartTime(false, false);

        // Update design variables
        optType_->update();

        // Restore time
        timeManip_.restoreTime();
    }
}


void Foam::unsteadyOptimisation::moveDesignVariables
(
    scalarField& direction
)
{
    if (shouldUpdateDesignVariables_)
    {
        // Updated grid should be written at the primal end time, to coinside
        // with the primal and adjoint solutions to be stored then.  The
        // (unknown) time index is not important at this point, so it remains
        // unchanged.
        timeManip_.storeTime();
        timeManip_.moveToAdjointStartTime(false, false);

        // Update design variables
        optType_->update(direction);

        // Restore time
        timeManip_.restoreTime();
    }
}


Foam::label Foam::unsteadyOptimisation::nActivePrimalSolvers()
{
    const dictionary& primalSolversDict =
        IOdictionary::subDict("primalSolvers");
    const wordList& primalSolverNames = primalSolversDict.toc();
    label n(0);
    for (const word& name : primalSolverNames)
    {
        if (primalSolversDict.subDict(name).getOrDefault<bool>("active", true))
        {
            n++;
        }
    }
    return n;
}


Foam::label Foam::unsteadyOptimisation::nActiveAdjointSolvers()
{
    const dictionary& adjointManagersDict =
        IOdictionary::subDict("adjointManagers");
    const wordList& adjointManagerNames = adjointManagersDict.toc();
    label n(0);
    for (const word& name : adjointManagerNames)
    {
        const dictionary& dict = adjointManagersDict.subDict(name);
        n += adjointSolverManager::nActiveAdjointSolvers(dict);
    }
    return n;
}


void Foam::unsteadyOptimisation::initialize()
{
    mesh_.moving(false); // Note: not sure if needed
    // Boolean defining if primal fields will be read from a custom time
    // provided by the user.
    Switch useCustomReadTime(true);
    Switch shouldWriteObjective(true);
    if (!localIOdictionary::empty())
    {
        // Possible scenarios for the solution of the primal/adjoint equations:
        // (a) one wishes to continue with the next optimization cycle, or
        // (b) has run the primal solver, has stored the flow fields
        // time-series and now wants to run adjoint,
        // (c) has already solved the adjoint equations and wants to integrate
        // them again based on the same primal fields and without proceeding to
        // the next optimization cycle, or
        // (d) wants to continue with the primal solver only.
        // For (a), nothing needs to be done here.
        // For (b), in the previous execution of the code only the flow
        // equations were integrated and, now, the adjoint equations are
        // about to be solved. Hence, isAdjointActive under uniform is false
        // and nActiveAdjointSolvers() > 0 and nActivePrimalSolvers() >= 0.
        // nActivePrimalSolvers() is allowed to take positive values since the
        // user may have stored the primal problem and now wants to run adjoint
        // only (nActivePrimalSolvers() = 0) or wants to run adjoint and, then,
        // continue with the optimization (nActivePrimalSolvers() > 0)
        // For (c), isAdjointActive under uniform is true,
        // nActiveAdjointSolvers() > 0 and nActivePrimalSolvers() = 0.
        // For (d), isAdjointActive under uniform is false,
        // nActiveAdjointSolvers() = 0 and nActivePrimalSolvers() > 0.
        // iOptCycles_ does not need  to be checked for the following reason.
        // Even if iOptCycles_ > 0 and the user wants to run the primal solver
        // only, then if
        // (i) the previous optimization cycle has just been completed case4 is
        // not activated and the manager proceeds as if a new optimization
        // cycle begins (a), whereas if
        // (ii) the primal solver has already been used, primal fields must be
        // initialized using the last available time-step

        label nPrimalSolvers = nActivePrimalSolvers();
        label nAdjointSolvers = nActiveAdjointSolvers();
        // Check if adjoint solvers were active in the previous execution
        Switch isAdjointActive =
            localIOdictionary::get<Switch>("isAdjointActive");

        Switch case2(!isAdjointActive && nAdjointSolvers > 0);
        Switch case3
        (
            isAdjointActive  && nAdjointSolvers > 0 && nPrimalSolvers == 0
        );
        Switch case4
        (
            !isAdjointActive && nAdjointSolvers == 0 && nPrimalSolvers > 0
        );
        DebugInfo
            << "case2 " << case2 << nl
            << "case3 " << case3 << nl
            << "case4 " << case4 << endl;
        if (case2 || case3)
        {
            // Start time of the primal is this time - the old span
            timeManip_.setStartTime
            (
                mesh_.time().value() - localIOdictionary::get<scalar>("span")
            );
            // Read startTimeIndex from dict
            timeManip_.setStartTimeIndex
                (localIOdictionary::get<label>("primalStartTimeIndex"));
            // Don't update design variables and avoid writing objectives
            // since we are running adjoint only
            shouldUpdateDesignVariables_ = false;
            shouldWriteObjective = false;
            // Don't solve the primal equations in the first go, even if they
            // are enabled
            solveFirstPrimalEqns_ = false;
        }
        else if (case4)
        {
            // Don't update design variables since we are running primal only
            shouldUpdateDesignVariables_ = false;
            // Disable editing readTime of the primal fields
            useCustomReadTime = false;
        }
        if (case3)
        {
            if (--iOptCycle_ < 0)
            {
                FatalErrorInFunction
                    << "iOptCycle = " << iOptCycle_ << ". Code should enter "
                    << "here only if both the primal and adjoint" << nl
                    << "equations have already been integrated for the" << nl
                    << "current optimization cycle and the user wants to" << nl
                    << "run the adjoint solver again" << endl
                    << exit(FatalError);
            }
        }
    }
    setSolvers(useCustomReadTime);
    // Deactivate writting objective functions in case that the primal
    // solution of the current optimization cycle has been stored and
    // the adjoint problem is about to be solved
    if (!shouldWriteObjective)
    {
        forAll(adjointSolverManagers_, amI)
        {
            adjointSolverManagers_[amI].setWrite(false);
            adjointSolverManagers_[amI].setWriteOption(IOobject::NO_WRITE);
        }
    }
    // Set averaging start time (if activated)
    forAll(primalSolvers_, psI)
    {
        primalSolvers_[psI].getVariablesSet().
            adjustAverageStartTime(timeManip_.startTime().value());
    }

//  if (debug)
//  {
        Info<< "startTime " << timeManip_.startTime() << nl
            << "startTimeIndex " << timeManip_.startTimeIndex() << nl
            << "primalEndTime " << timeManip_.primalEndTime() << nl
            << "span " << timeManip_.span() << endl;
//  }
}


void Foam::unsteadyOptimisation::setSolvers(Switch useCustomReadTime)
{
    // Read primal fields from the time prescribed in the optimisationDict,
    // Time index is not updated, but this should not be a problem
    dictionary& primalSolversDict =
        IOdictionary::subDict("primalSolvers");
    const wordList& primalSolverNames = primalSolversDict.toc();

    // Construct primal solvers
    primalSolvers_.setSize(primalSolverNames.size());
    forAll(primalSolvers_, solveri)
    {
        dictionary& solverDict =
            primalSolversDict.subDict(primalSolverNames[solveri]);
        if (primalSolvers_.size() > 1)
        {
            solverDict.add<bool>("useSolverNameForFields", true);
        }
        primalSolvers_.set
        (
            solveri,
            primalSolver::New
            (
                mesh_,
                managerType_,
                solverDict,
                useCustomReadTime
            )
        );
    }

    // Construct adjointSolverManagers
    const dictionary& adjointManagersDict =
        IOdictionary::subDict("adjointManagers");
    const wordList& adjointManagerNames = adjointManagersDict.toc();
    adjointSolverManagers_.setSize(adjointManagerNames.size());

    label nAdjointSolvers(0);
    label nActiveAdjointSolvers(0);
    bool overrideUseSolverName(adjointSolverManagers_.size() > 1);
    forAll(adjointSolverManagers_, manageri)
    {
        adjointSolverManagers_.set
        (
            manageri,
            new adjointSolverManager
            (
                mesh_,
                managerType_,
                adjointManagersDict.subDict(adjointManagerNames[manageri]),
                overrideUseSolverName
            )
        );
        nAdjointSolvers +=
            adjointSolverManagers_[manageri].nAdjointSolvers();
        nActiveAdjointSolvers +=
            adjointSolverManagers_[manageri].nActiveAdjointSolvers();
    }
    isAdjointActive_ = (nActiveAdjointSolvers > 0);

    // Sanity checks on the naming convention
    if (primalSolvers_.size() > 1)
    {
        for (const primalSolver& solveri : primalSolvers_)
        {
            if (!solveri.useSolverNameForFields())
            {
                FatalErrorInFunction
                    << "Multiple primal solvers are present but "
                    << "useSolverNameForFields is set to false in "
                    << "primal solver " << solveri.solverName() << nl
                    << "This is considered fatal."
                    << exit(FatalError);
            }
        }
    }

    if (nAdjointSolvers > 1)
    {
        for (const adjointSolverManager& amI : adjointSolverManagers_)
        {
            const PtrList<adjointSolver>& adjointSolvers = amI.adjointSolvers();
            for (const adjointSolver& asI : adjointSolvers)
            {
                if (!asI.useSolverNameForFields())
                {
                    FatalErrorInFunction
                        << "Multiple adjoint solvers are present but "
                        << "useSolverNameForFields is set to false in "
                        << "adjoint solver " << asI.solverName() << nl
                        << "This is considered fatal."
                        << exit(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unsteadyOptimisation::unsteadyOptimisation(fvMesh& mesh)
:
    optimisationManager(mesh),
    localIOdictionary
    (
        IOobject
        (
            "optimisation",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    iOptCycle_(localIOdictionary::getOrDefault<label>("optimisationCycle", 0)),
    nOptCycles_(mesh_.time().controlDict().get<label>("nOptimisationCycles")),
    isAdjointActive_(true),
    timeManip_
    (
        const_cast<unsteadyTimeManipulation&>
        (
            unsteadyTimeManipulation::New(mesh)
        )
    ),
    solveFirstPrimalEqns_(true)
{
    // Set primal and adjoint solvers
    initialize();
    optType_.reset
    (
        incompressible::optimisationType::New
        (
            mesh,
            IOdictionary::subDict("optimisation"),
            adjointSolverManagers_
        )
    );
    // Update source ptrs in all solvers to look at the source held in optType
    // Possible problem if mesh is adapted
    updateOptTypeSource();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::optimisationManager& Foam::unsteadyOptimisation::operator++()
{
    iOptCycle_++;

    if (!end())
    {
        Info<< "\n* * * * * * * * * * * * * * * * *" << endl;
        Info<< "Optimisation cycle " << iOptCycle_ << endl;
        Info<< "Time span ["
            << time_.value() << ","
            << time_.value() + timeManip_.span()
            << "]" << endl;
        Info<< "* * * * * * * * * * * * * * * * *\n" << endl;
    }
    timeManip_.newOptimisationCycle();

    return *this;
}


Foam::optimisationManager& Foam::unsteadyOptimisation::operator++(int)
{
    return operator++();
}


bool Foam::unsteadyOptimisation::checkEndOfLoopAndUpdate()
{
    if (update())
    {
        optType_->update();
    }
    return end();
}


bool Foam::unsteadyOptimisation::end()
{
    return iOptCycle_ > nOptCycles_;
}


bool Foam::unsteadyOptimisation::update()
{
    return (iOptCycle_ != 1 && !end());
}


void Foam::unsteadyOptimisation::solvePrimalEquations()
{
    // Solve all primal equations
    if (solveFirstPrimalEqns_)
    {
        optimisationManager::solvePrimalEquations();
    }
    solveFirstPrimalEqns_ = true;
}


void Foam::unsteadyOptimisation::solveAdjointEquations()
{
    optimisationManager::solveAdjointEquations();
    // Re-activate Switch to update the design variables at the
    // beggining of the next optimization cycle, if it has been
    // deactivated
    shouldUpdateDesignVariables_ = true;
    // Allow writing objectives at the next optimization cycle,
    // if deactivated
    forAll(adjointSolverManagers_, amI)
    {
        adjointSolverManagers_[amI].setWrite(true);
        adjointSolverManagers_[amI].setWriteOption(IOobject::AUTO_WRITE);
    }
}


bool Foam::unsteadyOptimisation::writeData(Ostream& os) const
{
    DebugInfo
        << "Calling unsteadyOptimisation::writeData" << endl;

    // If adjoint hasn't been solved, don't treat this cycle as complete
    label cycle(isAdjointActive_ ? iOptCycle_ : iOptCycle_ - 1);

    const label startTimeIndex = timeManip_.startTimeIndex();
    os.writeEntry("optimisationCycle", cycle);
    os.writeEntry("isAdjointActive", isAdjointActive_);
    os.writeEntry("primalStartTimeIndex", startTimeIndex);

    DebugInfo
        << "Writing data to localIOdictionary for continuation" << nl
        << tab << "optimisationCycle    " << cycle              << nl
        << tab << "isAdjointActive      " << isAdjointActive_   << nl
        << tab << "primalStartTimeIndex " << startTimeIndex     << endl;

    return true;
}


// ************************************************************************* //
