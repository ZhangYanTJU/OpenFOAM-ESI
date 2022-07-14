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
    storageParameters

Description
    Reads from dictionary and stores all parameters related to the storage
    of the primal fields

SourceFiles
    storageParameters.C

\*---------------------------------------------------------------------------*/

#ifndef storageParameters_H
#define storageParameters_H

#include "fvMesh.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class storageParameters Declaration
\*---------------------------------------------------------------------------*/

class storageParameters
{
    static label counter;

protected:


    // Protected Data

        const fvMesh& mesh_;

        const dictionary dict_;

        const wordList& allocatedFieldNames_;


        // Entries read from dictionary

            //- Entry defining if the simulation has a constant time-step or
            //- not. Added here to define whether the time-step must be updated
            //- when retrieving the variablesSet in checkPointing
            Switch adjustTimeStep_;

            //- Entry defining if all boundaries will be stored, regardless if
            //- they contain values that can be reproduced upon reconstruction.
            //  Defaults to false.
            Switch storeAllBoundaries_;

            //- Switch defining if the uniform boundaries must be compressed or
            //- not. If false, a single value (scalar or vector) will be kept
            //- for each of these boundaries. If true, these boundaries will be
            //- stored using the same algorithm with the internalField. If
            //- 'storeAllBoundaries' is set to true, this switch should have no
            //- effect
            Switch storeUniformBoundaries_;

            //- Algorithm to be used for the compression (ZFP, iPGD, etc)
            word algorithm_;

            //- Entry defining the algorithm variant to be used for the
            //- compression of each field.  Valid entries are the derived
            //- classes of the compressedGeometricField class.  The method is
            //- defined here once for all the objects to be created during the
            //- optimization (of course the method may differ for each primal
            //- field).
            wordList compressionMethod_;

            //- scalarList containing the times on which the primal fields must
            //- be written. Both before compression and after decompression.
            scalarList writeFieldTimes_;

            //- boolean controling if all unnecessary actions are avoided to
            //- benchmark the code
            bool timing_;

            //- boolean activated only if timing_ is set to false. If
            //- activated, the writeFieldTimes is ignored and the 'exact' and
            //- the reconstructed fields are written at all time-steps
            bool writeAll_;


        //- Small scalar value to be added in the denominator to avoid division
        //- by zero
        scalar tiny_;

        //- Return the vector of solved-for directions in mesh.
        //  1 indicates valid direction and 0 an invalid direction.
        const Vector<label> solDirs_;

        //- List containing the IDs of the first cell centers next to wall
        //- boundaries
        labelList wallList_;


    // Protected Member Functions

        //- Fills in wallList_
        void findCellsNextToWall();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        storageParameters(const storageParameters&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const storageParameters&) = delete;


public:

    TypeName("storageParameters");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            storageParameters,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& storageDict,
                const wordList& allocatedFieldNames
            ),
            (mesh, storageDict, allocatedFieldNames)
        );


    // Constructors

        //- Construct from components
        storageParameters
        (
            const fvMesh& mesh,
            const dictionary& storageDict,
            const wordList& allocatedFieldNames
        );


    // Selectors

        //- Return a reference to the selected storageParameters
        static autoPtr<storageParameters> New
        (
            const fvMesh& mesh,
            const dictionary& storageDict,
            const wordList& allocatedFieldNames
        );


    //- Destructor
    virtual ~storageParameters() = default;


    // Member functions

        virtual void initialize();

        virtual void reset();

        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        inline const Switch& adjustTimeStep() const
        {
            return adjustTimeStep_;
        }

        inline const Switch& storeAllBoundaries() const
        {
            return storeAllBoundaries_;
        }

        inline const Switch& storeUniformBoundaries() const
        {
            return storeUniformBoundaries_;
        }

        inline const word& algorithm() const
        {
            return algorithm_;
        }

        inline const dictionary& dict() const
        {
            return dict_;
        }

        inline const scalarList& writeFieldTimes() const
        {
            return writeFieldTimes_;
        }

        inline bool timing() const
        {
            return timing_;
        }

        inline bool writeAll() const
        {
            return writeAll_;
        }

        inline const Vector<label>& solDirs() const
        {
            return solDirs_;
        }

        inline const wordList& compressionMethod() const
        {
            return compressionMethod_;
        }

        inline scalar tiny() const
        {
            return tiny_;
        }

        inline const labelList& wallList() const
        {
            return wallList_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
