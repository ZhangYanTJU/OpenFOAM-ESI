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


Class
    Foam::pimple

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef pimple_H
#define pimple_H

#include "incompressiblePrimalSolver.H"
#include "PIMPLEControl.H"
#include "IOMRFZoneList.H"
#include "objective.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class pimple Declaration
\*---------------------------------------------------------------------------*/

class pimple
:
    public incompressiblePrimalSolver
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pimple(const pimple&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const pimple&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<PIMPLEControl> solverControl_;

        //- Boolean defining if primal fields will be initialized from
        //  the current or a costum time-step provided by the user
        Switch shouldResetMeanFields_;

        //- Reference to incompressibleVars
        //  Used for convenience and to avoid repetitive dynamic_casts
        //  Same as getIncoVars()
        incompressibleVars& incoVars_;

        //- MRF zones
        IOMRFZoneList MRF_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

        //- List of objectives related to this primal solver
        List<objective*> objectives_;

        //- Boolean and scalar entries for readDyMControls and
        //  setDeltaT functions
        bool adjustTimeStep_;
        scalar maxCo_;
        scalar maxDeltaT_;
        scalar CoNum_;
        /*bool correctPhi_;
        bool checkMeshCourantNo_;
        bool moveMeshOuterCorrectors_;*/

        //- velocity field Uf (set only if required)
        autoPtr<surfaceVectorField> Uf_;

        scalar deltaT0_;

        //- Time-index at the beggining of the primal solution
        label primalStartTimeIndex_;

        //- List with deltaT values during the solution of the primal problem
        DynamicList<scalar> deltaTList_;

        //- Time manipulation, common for all primal and adjoint solvers
        unsteadyTimeManipulation& timeManip_;


    // Protected Member Functions

        //- Allocate incompressibleVars and return reference to be used for
        //- convenience in the rest of the class.
        incompressibleVars& allocateVars(scalar readTime);

        //- Read only the control parameters used by setDeltaT
        void readTimeControls();

        //- Compute Courant No
        void CourantNo();

        //- Compute continuity errors
        void continuityErrors();

        //- Reset time-step to maintain a constant maximum courant Number.
        //  Reduction of time-step is immediate, but increase is damped to
        //  avoid unstable oscillations.
        void setDeltaT();

        //- Set the initial timestep corresponding to the timestep
        //  adjustment algorithm in setDeltaT but only if it would
        //  reduce the timestep.
        void setInitialDeltaT();

        //- Increase time
        void setTime();

        //- Store the current deltaT
        void storeDeltaT();

        //- Scratch deltaTList_ and deltaT0_
        void scratchDeltaTHistory();

        /*
        //- Read among other the control parameters used by setDeltaT
        void readDyMControls();

        //- Calculates and outputs the mean and maximum Courant Numbers.
        void meshCourantNo();

        //- Flux correction functions to ensure continuity
        void correctPhi();

        //- Creates and initialises the velocity field Uf if required.
        void createUfIfPresent();
        */


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pimple");


    // Constructors

        //- Construct from mesh and dictionary
        pimple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            Switch useCustomReadTime = false
        );


    //- Destructor
    virtual ~pimple() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict) override;

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter() override;

            //- Main control loop
            virtual void solve() override;

            //- Looper (advances iters, time step)
            virtual bool run();
            virtual bool loop() override;

            //- Restore initial field values if necessary
            virtual void restoreInitValues() override;

            //- Functions to be called before loop
            virtual void preLoop() override;

            //- Functions to be called after loop
            virtual void postLoop() override;

            //- Functions to be called after lineSearch has converged
            virtual void postLineSearch() override;

            //- Write average iteration
            virtual bool writeData(Ostream& os) const override;


        // Info to the AdjointSolver

            //- Return a scalarList with all deltaT of the primal solution
            const DynamicList<scalar>& getDeltaTList() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
