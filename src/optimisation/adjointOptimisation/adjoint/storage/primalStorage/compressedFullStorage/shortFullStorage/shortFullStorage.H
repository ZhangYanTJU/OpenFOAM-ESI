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
    Foam::shortFullStorage

Description
    Stores all instances of the primal flow fields.
    In constrast to full storage, parts of the fields that can be recomputed
    (fixedValue boundary conditions, processor values etc) are not stored,
    reducing slightly the memory footprint

SourceFiles
    shortFullStorage.C

\*---------------------------------------------------------------------------*/

#ifndef shortFullStorage_H
#define shortFullStorage_H

#include "compressedFullStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class shortFullStorage Declaration
\*---------------------------------------------------------------------------*/

class shortFullStorage
:
    public compressedFullStorage
{
protected:

    // Protected Data

        label& iPtr_;


    // Protected Member Functions

        //- Gather metrics from all processors if timing == true.
        //  If timing == false, the reduce operation has already been deployed.
        void gatherMetrics();

        //- Set file to write metrics - one line per optimisation cycle
        void setMetricsFilesPtr();

        inline void writeToFile();

        //- Function to calculate the compression metrics of the whole primal
        //- solution
        void calcStorageMetrics();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        shortFullStorage(const shortFullStorage&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const shortFullStorage&) = delete;


public:

    TypeName("shortFullStorage");

    // Constructors

        //- Construct from components
        shortFullStorage
        (
            primalSolver& primalSolverObj,
            const dictionary storageDict
        );


    //- Destructor
    virtual ~shortFullStorage() = default;


    // Member functions

        //- Free memory used to store the primal time-series
        virtual void scratch() override;

        //- Store current copy of the variablesSet (could be lossy)
        virtual void storeVariables() override;

        //- Retrieve primal-based quantities into the active variablesSet
        virtual void retrieveVariables() override;

        //- Function to write the storage metrics of the whole primal solution
        void storageMetrics() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
