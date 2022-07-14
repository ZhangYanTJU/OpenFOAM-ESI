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
    Foam::binomialCheckPointing

Description
    Class implementing the binomial check-pointing technique for the unsteady
    adjoint optimisation

SourceFiles
    binomialCheckPointing.C

\*---------------------------------------------------------------------------*/

#ifndef binomialCheckPointing_H
#define binomialCheckPointing_H

#include "primalStorage.H"
#include "storageParameters.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "checkPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class binomialCheckPointing Declaration
\*---------------------------------------------------------------------------*/

class binomialCheckPointing
:
    public primalStorage
{

protected:

    // Protected Data

        //- Counter of the total number of time-steps performed to
        //- complete the optimization cycle
        //  primal + adjoint + primal recomputations due to checkPointing
        label totalSteps_;

        //- Counter of the number of recomputation time-steps
        //- performed to complete the optimization cycle
        label recomputationSteps_;

        //- Storage parameters
        autoPtr<storageParameters> storageParams_;

        //- Reseting deltaT causes slight disrepancies in comparison to the
        //- solution of shortFullStorage.
        //  Reseting applied only with non-constant time-step. Entry taken from
        //  storageParameters to avoid the possibility of the user changing
        //  runTime the adjustTimeStep of controlDict
        const Switch& adjustTimeStep_;

        //- Number of chechPoints to be stored
        label nCheckPoints_;

        //- Accumulated storage metrics
        scalarList storageMetrics_;

        //- All stored checkPoints
        PtrList<checkPoint> checkPoints_;

        //- Position of each checkPoint, as stored in chronological order,
        //- in the checkPoints list
        List<label> indirectAddressing_;

        dimensionedScalar startingTime_;

        label activeCheckPoints_;

        label& iPtr_;

        mutable autoPtr<OFstream> metricsFilePtr_;

        mutable autoPtr<OFstream> checkPointsFilePtr_;


    // Protected Member Functions

        //- Gather metrics from all processors if timing == true.
        //  If timing == false, the reduce operation has already been deployed.
        void gatherMetrics();

        //- Set file to write metrics - one line for each optimisation cycle
        void setMetricsFilesPtr();

        //- Set file to write checkPoints
        void setCheckPointsFilesPtr();

        inline void writeToFile();

        //- Calculates the compression metrics of the whole primal solution
        virtual void calcStorageMetrics();

        //- Start process for each host
        labelList hostStartProcess(wordList& hostNames) const;

        //- Get checkpoint size. The size of the first checkpoint size is not
        //- dependable, since it may include only the fields of the current
        //- time-step and not the oldTimes.  Instead, we solve for one
        //- time-step and use its size as a reference
        scalar checkPointSize();

        //- Print checkpoint size and estimate number of checkpoints that
        //- can be allocated based on the available memory at run time
        unsigned long long maxNumberOfCheckPoints();

        void initialize();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        binomialCheckPointing(const binomialCheckPointing&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const binomialCheckPointing&) = delete;


public:

    TypeName("binomialCheckPointing");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binomialCheckPointing,
            dictionary,
            (
                primalSolver& primalSolverObj,
                const dictionary storageDict
            ),
            (primalSolverObj, storageDict)
        );


    // Constructors

        //- Construct from components
        binomialCheckPointing
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    // Selectors

        //- Allocate new binomialCheckPointing
        static autoPtr<binomialCheckPointing> New
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    //- Destructor
    virtual ~binomialCheckPointing() = default;


    // Member functions

        //- Release all checkpoints
        virtual void scratch() override;

        //- Store current variable set as a checkpoint if needed
        virtual void storeVariables() override;

        virtual void storeInitialVariables() override;

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
