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
    Foam::fullStorage

Description
    Stores all instances of the primal flow fields.
    Can lead to huge memory requirements if used with large meshes and many
    time-steps

SourceFiles
    fullStorage.C

\*---------------------------------------------------------------------------*/

#ifndef fullStorage_H
#define fullStorage_H

#include "primalStorage.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fullStorage Declaration
\*---------------------------------------------------------------------------*/

class fullStorage
:
    public primalStorage
{
protected:

    // Protected Data

        PtrList<variablesSet> allPrimalVars_;

        label& iPtr_;


    // Protected Member Functions

        //- Start process for each host
        labelList hostStartProcess(wordList& hostNames) const;

        //- Print time instance size and estimate number of instances that
        //- can be allocated based on the available memory at run time
        label maxNumberOfInstances() const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fullStorage(const fullStorage&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fullStorage&) = delete;


public:

    TypeName("fullStorage");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fullStorage,
            dictionary,
            (
                primalSolver& primalSolverObj,
                const dictionary storageDict
            ),
            (primalSolverObj, storageDict)
        );


    // Constructors

        //- Construct from components
        fullStorage
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<fullStorage> New
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    //- Destructor
    virtual ~fullStorage() = default;


    // Member functions

        virtual void scratch() override;

        virtual void storeVariables() override;

        virtual void retrieveVariables() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
