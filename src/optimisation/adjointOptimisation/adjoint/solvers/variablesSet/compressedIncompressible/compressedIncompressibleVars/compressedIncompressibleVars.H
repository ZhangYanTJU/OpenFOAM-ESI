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
    Foam::compressedIncompressibleVars

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef compressedIncompressibleVars_H
#define compressedIncompressibleVars_H

#include "incompressibleVars.H"
#include "compressedGeometricFields.H"
#include "storageParameters.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class compressedIncompressibleVars Declaration
\*---------------------------------------------------------------------------*/

class compressedIncompressibleVars
:
    public variablesSet
{
protected:

    // Protected data

        //- Vars seen by the solver
        incompressibleVars& incoVars_;

        //- Solver control related to the base variables
        const solverControl& solverControl_;

        storageParameters& storageParams_;

        //- Current S/N of the time-step that is compressed in the present
        //- object.
        //  Useful when more than one time-steps are compressed
        //  in the same compressedIncompressibleVars (e.g. iPGD algorithm)
        label k_;

        label kCopy_;

        //- Compressed fields including turbulence
        autoPtr<compressedVolScalarField> p_;
        autoPtr<compressedSurfaceScalarField> phi_;
        autoPtr<compressedVolVectorField> U_;
        PtrList<compressedVolScalarField> RASModelVars_;

        //- Time value and index
        label timeIndex_;

        scalar timeValue_;

        //- List that holds the compression metrics (initial size, compressed
        //- size if compression is applied, etc) for all fields of the current
        //- time step.
        scalarList storageMetrics_;

        //- Variables used for writing the compression metrics
        fileName storageFolder_;
        wordList names_;
        mutable PtrList<OFstream> storageFilePtr_;


    // Protected Member Functions

        //- Set compressed fields including turbulence
        void setFields();

        //- Calculate and write compression metrics for a
        //- compressedGeometricField
        template<class Type, template<class> class PatchField, class GeoMesh>
        void calculateAndWrite
        (
            const compressedGeometricField<Type, PatchField, GeoMesh>& compField,
            const word& name,
            label& i
        );

        //- Function to write the compressed size, the compression ratio, etc.
        virtual void writeLog(const word&) const;

        //- Function to write compression metrics to file
        virtual void write(const label);

        //- Set folder
        void makeFolder();

        //- Set File to store compression metrics
        virtual void setStorageFilesPtr();


public:


    // Static Data Members

        //- Run-time type information
        TypeName("compressedIncompressibleVars");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressedIncompressibleVars,
            dictionary,
            (
                incompressibleVars& vs,
                storageParameters& storageParams
            ),
            (vs, storageParams)
        );


    // Constructors

        compressedIncompressibleVars
        (
            incompressibleVars& vs,
            storageParameters& storageParams
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<compressedIncompressibleVars> New
        (
            incompressibleVars& vs,
            storageParameters& storageParams
        );


    //- Destructor
    virtual ~compressedIncompressibleVars() = default;


    // Member Functions

        inline label timeIndex() const
        {
            return timeIndex_;
        }

        inline label& timeIndex()
        {
            return timeIndex_;
        }

        inline scalar timeValue() const
        {
            return timeValue_;
        }

        inline scalar& timeValue()
        {
            return timeValue_;
        }

        //- Function to calculate and write the compression metrics
        virtual void calculateStorageMetrics();

        virtual void compress();

        virtual void decompress
        (
            incompressibleVars& vars
        );

        void decompress();

        //- Function to return the compression metrics (initial size,
        //- compressed size if compression is applied, etc)
        virtual const scalarList& storageMetrics() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "compressedIncompressibleVarsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
