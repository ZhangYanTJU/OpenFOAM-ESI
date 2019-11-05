/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::fa::externalFileSource

Group
    grpFvOptionsSources

Description


Usage

    \verbatim
    externalFileSourceCoeffs
    {
        pressure
        {
            type                                externalFileSource;
            fieldName                           ws_ceilingShell;
            tableName                           pressude.dat
            active                              true;
            region                              ceilingShell;
            selectionMode                       all;
        }
    }
    \endverbatim


    Where:
    \table
        Property     | Description               | Required  | Default value
        fieldName    | Name of the field         | yes
        tableName    | Name of the data table    | yes
        region       | 2D Region to be applied   | yes
    \endtable


SourceFiles
    externalFileSource.C

SeeAlso
    Foam::PatchFunction1Types

\*---------------------------------------------------------------------------*/

#ifndef fa_externalFileSource_H
#define fa_externalFileSource_H

#include "faOption.H"
#include "areaFields.H"
#include "faceSetOption.H"
#include "MappedFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class externalFileSource Declaration
\*---------------------------------------------------------------------------*/

class externalFileSource
:
    public faceSetOption
{
    // Private data

        //- Name of the field to apply this source
        word fieldName_;
        
        //- Name of the table
        word tableName_;
        
        //- Mapped data from file
        autoPtr<PatchFunction1Types::MappedFile<scalar>> value_;
 
        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- No copy construct
        externalFileSource(const externalFileSource&) = delete;

        //- No copy assignment
        void operator=(const externalFileSource&) = delete;

        //- Initialise the electrical conductivity field
        template<class Type>
        void initialiseSigma
        (
            const dictionary& dict,
            autoPtr<Function1<Type>>& sigmaVsTPtr
        );

        //- Update the electrical conductivity field
        template<class Type>
        const GeometricField<Type, faPatchField, areaMesh>&
        updateSigma(const autoPtr<Function1<Type>>& sigmaVsTPtr) const;


public:

    //- Runtime type information
    TypeName("externalFileSource");


    // Constructors

        //- Construct from explicit source name and mesh
        externalFileSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~externalFileSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
