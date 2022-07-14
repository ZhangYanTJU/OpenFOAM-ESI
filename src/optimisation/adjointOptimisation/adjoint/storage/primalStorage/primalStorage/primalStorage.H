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
    Foam::primalStorage

Description
    Abstract base class managing the storage of primal flow fields, to be later
    used for the solution of unsteady flow equations

SourceFiles
    primalStorage.C

\*---------------------------------------------------------------------------*/

#ifndef primalStorage_H
#define primalStorage_H

#include "Time.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "variablesSet.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

// Forward declaratrions
class primalSolver;

/*---------------------------------------------------------------------------*\
                        Class primalStorage Declaration
\*---------------------------------------------------------------------------*/

class primalStorage
{
    static label counter;

protected:

    // Protected Data

        //- Const reference to time
        const Time& time_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Reference to the primal solver used to create the time-series
        //- being stored
        primalSolver& primalSolver_;

        //- Input dict
        dictionary storageDict_;

        //- Variables set of the primal solver
        variablesSet& variablesSet_;

        //- List of useful pointers in the flow time-series
        //  For instance, time-steps holding checkpoints, etc
        labelList storagePtrs_;

        //- Copy of the pointer's list
        //  Usefull for reseting in many of many primal solvers
        autoPtr<labelList> storagePtrsCopy_;

        //- Output folder
        fileName storageFolder_;


    // Protected Member Functions

        //- Set Folder
        virtual void makeFolder();

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        primalStorage(const primalStorage&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const primalStorage&) = delete;


public:

    TypeName("primalStorage");

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            primalStorage,
            primalStorage,
            (
                primalSolver& primalSolverObj,
                const dictionary storageDict
            ),
            (primalSolverObj, storageDict)
        );


    // Constructors

        //- Construct from components
        primalStorage
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<primalStorage> New
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    // Destructor

        //- Virtual destructor
        virtual ~primalStorage() = default;


    // Member functions

        //- Get list with pointers to storage instances
        virtual const labelList& storagePtrs() const;

        //- Copy pointers to storage instances to a backup variable
        virtual void copyStoragePtrs();

        //- Rewind pointers to storage instance to the values set after the
        //- solution of the primal equations
        virtual void rewindStoragePtrs();

        //- Free memory used to store the primal time-series
        virtual void scratch() = 0;

        //- Store current copy of the variablesSet (could be lossy)
        virtual void storeVariables() = 0;

        virtual void storeInitialVariables();

        //- Retrieve primal-based quantities into the active variablesSet
        virtual void retrieveVariables() = 0;

        //- Retrive certain variables before the adjoint loop
        //  (e.g. from files, for restarted cases)
        virtual void preLoopRetrieveVariables();

        //- Operations performed before the beginning of the adjoint loop
        virtual void preAdjointLoop();

        //- Operations performed after the end of the adjoint loop
        virtual void postAdjointLoop();

        virtual void storageMetrics();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
