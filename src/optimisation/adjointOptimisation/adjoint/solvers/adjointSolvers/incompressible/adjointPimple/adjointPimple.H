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
    Foam::adjointPimple

Description
    Solution of the incompressible, unsteady adjoint equations using the PIMPLE
    algorithm

\*---------------------------------------------------------------------------*/

#ifndef adjointPimple_H
#define adjointPimple_H

#include "incompressibleAdjointSolver.H"
#include "PIMPLEControl.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "adjointSensitivityIncompressible.H"
#include "primalStorage.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class adjointPimple Declaration
\*---------------------------------------------------------------------------*/

class adjointPimple
:
    public incompressibleAdjointSolver
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        adjointPimple(const adjointPimple&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const adjointPimple&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<PIMPLEControl> solverControl_;

        //- Reference to incompressibleAdjointVars
        //  Used for convenience and to avoid repetitive dynamic_casts
        //  Same as getAdjointVars()
        incompressibleAdjointVars& adjointVars_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

        //- Sensitivity Derivatives engine
        autoPtr<incompressible::adjointSensitivity> adjointSensitivity_;

        //- Reference to primal storage mechanism
        autoPtr<primalStorage>& primalStorage_;

        //- labelList containing all the time-step dependent pointers necessary
        //  for the sucessfull retrieving of the variableSet
        labelList storagePtrs_;

        //- Time instances and indices at the beginning and end of each
        //- optimisation loop
        dimensionedScalar primalEndTime_;

        //- Reference to the deltaTList of the primal solver
        const DynamicList<scalar>& deltaTList_;

        //- Time manipulation, common for all primal and adjoint solvers
        unsteadyTimeManipulation& timeManip_;


    // Protected Member Functions

        //- Allocate incompressibleAdjointVars and return reference to be used
        //- for convenience in the rest of the class.
        incompressibleAdjointVars& allocateVars();

        //- Reset time-step to maintain a constant maximum courant Number.
        //  Reduction of time-step is immediate, but increase is damped to
        //  avoid unstable oscillations.
        void setDeltaT();

        //- Decrease time
        void setTime();

        //- Compute Courant No
        void CourantNo() const;

        //- Compute continuity errors
        void continuityErrors();

        //- Accumulate sensitivities from the various time instances
        void accumulateSensitivities();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointPimple");


    // Constructors

        //- Construct from mesh and dictionary
        adjointPimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~adjointPimple() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool run();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Functions to be called before loop
            virtual void preLoop();

            //- Functions to be called after loop
            virtual void postLoop();

            //- Functions to be called after the convergence of linearSearch
            virtual void postLineSearch();

            //- Compute sensitivities of the underlaying objectives
            virtual void computeObjectiveSensitivities();

            //- Grab a reference to the computed sensitivities
            virtual const scalarField& getObjectiveSensitivities();

            //- Clears the sensitivity field known by the adjoint solver
            //- and zeros sensitivities constituents
            virtual void clearSensitivities();

            //- Return the base sensitivity object
            virtual sensitivity& getSensitivityBase();

            //- Update primal based quantities
            //- related to the objective functions.
            //  Also writes the objective function values to files.
            //  Written here and not in the postLoop function of the primal
            //  to make sure we don't pollute the objective files with
            //  objectives of non-converged linearSearch iterations
            virtual void updatePrimalBasedQuantities();

            //- Write average iteration
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
