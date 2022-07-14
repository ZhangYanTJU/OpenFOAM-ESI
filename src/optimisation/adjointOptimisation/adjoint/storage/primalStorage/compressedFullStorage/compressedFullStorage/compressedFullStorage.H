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
    Foam::compressedFullStorage

Description
    Compresses and stores all instances of the primal flow fields

SourceFiles
    compressedFullStorage.C

\*---------------------------------------------------------------------------*/

#ifndef compressedFullStorage_H
#define compressedFullStorage_H

#include "primalStorage.H"
#include "compressedIncompressibleVars.H"
#include "storageParameters.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressedFullStorage Declaration
\*---------------------------------------------------------------------------*/

class compressedFullStorage
:
    public primalStorage
{
protected:

    // Protected Data

        PtrList<compressedIncompressibleVars> allPrimalVars_;

        scalarList storageMetrics_;

        //- Counter of the total number of time-steps performed to
        //- complete the optimization cycle.
        //  Equal to primal + adjoint + primal recomputations
        //  in case of checkPointing strategies
        label totalSteps_;

        //- Counter of the number of recomputation time-steps
        //- performed to complete the optimization cycle
        label recomputationSteps_;

        //- Storage parameters
        autoPtr<storageParameters> storageParams_;

        mutable autoPtr<OFstream> metricsFilePtr_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        compressedFullStorage(const compressedFullStorage&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const compressedFullStorage&) = delete;


public:

    TypeName("compressedFullStorage");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressedFullStorage,
            dictionary,
            (
                primalSolver& primalSolverObj,
                const dictionary storageDict
            ),
            (primalSolverObj, storageDict)
        );


    // Constructors

        //- Construct from components
        compressedFullStorage
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<compressedFullStorage> New
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    //- Destructor
    virtual ~compressedFullStorage() = default;


    // Member functions

        //- Free memory used to store the primal time-series
        virtual void scratch() override;

        //- Store current copy of the variablesSet (could be lossy)
        virtual void storeVariables() override;

        //- Retrieve primal-based quantities into the active variablesSet
        virtual void retrieveVariables() override;

        //- Operations performed after the end of the adjoint loop
        virtual void postAdjointLoop() override;

        //- Function to write the storage metrics of the whole primal solution
        virtual void storageMetrics() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
