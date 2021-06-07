/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Namespace
    Foam::binModels

Description
    A namespace for various bin model implementations.

Class
    Foam::binModel

Description
    Abstract base class for bin models to handle bin
    characteristics for \c force and \c forceCoeffs function objects.

See also
  - Foam::binModel
  - Foam::binModels::noBin
  - Foam::binModels::singleDirectionUniformBin
  - Foam::binModels::uniformBin

SourceFiles
    binModel.C
    binModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef binModel_H
#define binModel_H

#include "dictionary.H"
#include "HashSet.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                            Class binModel Declaration
\*---------------------------------------------------------------------------*/

class binModel
{
protected:

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flag to include porosity effects
        bool porosity_;

        //- Total number of bins (internal)
        label nBin_;

        //- Names of operand patches
        labelHashSet patchSet_;


        // File streams

            //- File stream for forces in bins
            autoPtr<OFstream> forceBinFilePtr_;

            //- File stream for moments in bins
            autoPtr<OFstream> momentBinFilePtr_;

            //- File streams for force and moment coefficients in bins
            PtrList<OFstream> coeffBinFilePtrs_;


    // Generated Methods

        //- No copy construct
        binModel(const binModel&) = delete;

        //- No copy assignment
        void operator=(const binModel&) = delete;


public:

    //- Runtime type information
    TypeName("binModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binModel,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, dict, mesh)
        );


    // Selectors

        //- Return a reference to the selected bin model
        static autoPtr<binModel> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        binModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~binModel() = default;


    // Member Functions

        // Access

            //- Read the dictionary
            virtual bool read(const dictionary& dict);

            //- Return the total number of bins
            label nBin() const noexcept
            {
                return nBin_;
            };


        // Evaluate

            //- Initialise bin properties
            virtual void initialise() = 0;

            //- Apply bins to force and moment containers
            virtual void applyBins
            (
                List<Field<vector>>& force,
                List<Field<vector>>& moment,
                const vectorField& d,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT
            ) = 0;

            //- Apply bins to force and moment
            //- containers, and include porosity effects
            virtual void applyBins
            (
                List<Field<vector>>& force,
                List<Field<vector>>& moment,
                const vectorField& d,
                const vectorField& Md,
                const vectorField& fN,
                const vectorField& fT,
                const vectorField& fP
            ) = 0;


        // Write

            //- Create binned-data files
            virtual void createBinnedDataFiles() = 0;

            //- Create a binned-data file - coefficients
            virtual void createBinnedDataFile
            (
                const word& title,
                const label id
            ) = 0;

            //- Write binned force and moment data to files
            virtual void writeBinnedDataFiles
            (
                const List<Field<vector>>& force,
                const List<Field<vector>>& moment,
                const coordSystem::cartesian& coordSys
            ) = 0;

            //- Write binned force and moment coefficient data to a file
            virtual void writeBinnedDataFile
            (
                List<Field<scalar>> coeff,
                const label id
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
