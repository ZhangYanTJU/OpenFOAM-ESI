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
    Foam::piso

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef piso_H
#define piso_H

#include "incompressiblePrimalSolver.H"
#include "PISOControl.H"
#include "IOMRFZoneList.H"
#include "objective.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class piso Declaration
\*---------------------------------------------------------------------------*/

class piso
:
    public incompressiblePrimalSolver
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        piso(const piso&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const piso&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<PISOControl> solverControl_;

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

        //- Time manipulation, common for all primal and adjoint solvers
        unsteadyTimeManipulation& timeManip_;


    //- Protected Member Functions

        //- Print flow fields if in debug mode
        void printFields();

        //- Allocate incompressibleVars and return reference to be used for
        //- convenience in the rest of the class.
        incompressibleVars& allocateVars(scalar readTime);

        //- Compute Courant No
        void CourantNo() const;

        //- Compute continuity errors
        void continuityErrors();

public:

    // Static Data Members

        //- Run-time type information
        TypeName("piso");


    // Constructors

        //- Construct from mesh and dictionary
        piso
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            Switch useCustomReadTime = false
        );


    //- Destructor
    virtual ~piso() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter() override;

            //- Main control loop
            virtual void solve() override;

            //- Looper (advances iters, time step)
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
