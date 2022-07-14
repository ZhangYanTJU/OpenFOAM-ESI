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


Class
    Foam::unsteadyOptimisation

Description
    Iterate the optimisation cycles. For steady state opt, this coinsides
    with evolving Time

SourceFiles
    unsteadyOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef unsteadyOptimisation_H
#define unsteadyOptimisation_H

#include "optimisationManager.H"
#include "localIOdictionary.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class unsteadyOptimisation Declaration
\*---------------------------------------------------------------------------*/

class unsteadyOptimisation
:
    public optimisationManager,
    public localIOdictionary
{
private:

    // Private Member Data

        //- Current optimisation cycle
        label iOptCycle_;

        //- Number of optimisation cycles
        label nOptCycles_;

        //- Switch defining if at least one adjoint solver is active
        Switch isAdjointActive_;

        //- Time manipulation, common for all primal and adjoint solvers
        unsteadyTimeManipulation& timeManip_;

        //- Boolean to treat a special case
        //  If adjoint was not solved in the previous code execution, but it is
        //  solved now, together with the primal equations, avoid solving the
        //  primal equations the first time
        bool solveFirstPrimalEqns_;


    // Private Member Functions

        //- Update optimisationType source for all primal and adjoint solvers
        void updateOptTypeSource();

        //- Disallow default bitwise copy construct
        unsteadyOptimisation(const unsteadyOptimisation&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const unsteadyOptimisation&) = delete;


protected:

    // Protected Member Functions

        //- Reset time to primal time after a failed line search
        virtual void resetTime() override;

        virtual void moveDesignVariables() override;

        virtual void moveDesignVariables
        (
            scalarField& direction
        ) override;

        //- Return number of active primal and adjoint solvers without
        //- constructing any of the solvers.
        //  It is calculated through reading dictionary entries. Warning: It is
        //  designed to be used by initialize() function only. In this case no
        //  mismatch can occur if the optimisationDict is changed during
        //  runTime, since the primal or adjoint solvers are built directly
        //  after calling this function. If the number of active adjoint
        //  solvers is required elsewhere, use the nActiveAdjointSolvers()
        //  function of each adjointSolverManager.
        label nActivePrimalSolvers();

        label nActiveAdjointSolvers();

        virtual void initialize() override;

        void setSolvers(Switch editPrimalReadTime);

public:

    //- Runtime type information
    TypeName("unsteadyOptimisation");


    // Constructors

        //- Construct from components
        unsteadyOptimisation(fvMesh& mesh);


    // Destructor
    virtual ~unsteadyOptimisation() = default;


    // Member Functions

        //- Prefix increment
        virtual optimisationManager& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual optimisationManager& operator++(int);

        //- Return true if end of optimisation run
        //  Also, updates the design variables if needed
        virtual bool checkEndOfLoopAndUpdate();

        //- Return true if end of optimisation run
        virtual bool end();

        //- Whether to update the design variables
        virtual bool update();

        //- Solve all primal equations
        virtual void solvePrimalEquations() override;

        //- Solve all adjoint equations
        virtual void solveAdjointEquations() override;

        //- write necessary data for continuation
        virtual bool writeData(Ostream&) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
