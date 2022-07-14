/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::incompressibleVars

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef incompressibleVars_H
#define incompressibleVars_H

#include "variablesSet.H"
#include "fvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "RASModelVariables.H"
#include "solverControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleVars Declaration
\*---------------------------------------------------------------------------*/

class incompressibleVars
:
    public variablesSet
{
protected:

    // Protected data

        //- Reference to the solverControl of the solver allocating the fields
        solverControl& solverControl_;

        //- Fields involved in the solution of the incompressible NS equations
        autoPtr<volScalarField> pPtr_;
        autoPtr<volVectorField> UPtr_;
        autoPtr<surfaceScalarField> phiPtr_;
        autoPtr<singlePhaseTransportModel> laminarTransportPtr_;
        autoPtr<incompressible::turbulenceModel> turbulence_;
        autoPtr<incompressible::RASModelVariables> RASModelVariables_;
        //autoPtr<volScalarField> TPtr_;

        //- Keep a copy of the initial values if necessary.
        //  If oldTimes are present, keep them in separate GeometricFields to
        //  avoid been overwritten upon calling GeometricField::storeOldTimes
        PtrList<volScalarField> pInitPtr_;
        PtrList<volVectorField> UInitPtr_;
        PtrList<surfaceScalarField> phiInitPtr_;

        //- Manage mean fields. Turbulent mean fields are managed
        //- in RASModelVariables
        autoPtr<volScalarField> pMeanPtr_;
        autoPtr<volVectorField> UMeanPtr_;
        autoPtr<surfaceScalarField>phiMeanPtr_;

        //- Update boundary conditions upon construction.
        //  Useful for cases in which information has been lost on boundary
        //  (e.g. averaging, redistribution)
        bool correctBoundaryConditions_;

        //- Names of the primal fields that are allocated.
        //  The size of the list is properly adjusted
        wordList allocatedFieldNames_;

        //- Return the vector of solved-for directions in mesh.
        //  1 indicates valid direction and 0 an invalid direction.
        const Vector<label> solDirs_;

        //- Should the solver fields be written
        bool writeFields_;


    // Protected Member Functions

        //- Read fields and set turbulence
        void setFields();

        //- Set initial fields if necessary
        void setInitFields();

        //- Set mean fields if necessary
        void setMeanFields();

        //- Rename turbulence fields if necessary
        void renameTurbulenceFields();

        //- Update boundary conditions of mean-flow
        void correctNonTurbulentBoundaryConditions();

        //- Update boundary conditions of turbulent fields
        void correctTurbulentBoundaryConditions();

        //- No copy assignment
        void operator=(const incompressibleVars&);

        //- Function to fill the wordList containing the names of the primal
        //- fields that are allocated.
        void storeAllocatedFieldNames();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("incompressibleVars");


    // Constructors

        //- Construct from mesh
        incompressibleVars
        (
            fvMesh& mesh,
            solverControl& SolverControl,
            scalar readTime = -1
        );

        //- Copy constructor
        //  turbulence_ and laminarTransport_ are uninitialized since there is
        //  no clear way (and possibly no usage) of cloning them.
        //  The resulting incompressibleVars is hence incomplete.
        //  Additionally, copied fields have the original name, appended with
        //  the timeName, to avoid database collisions.
        //  To be used as storing space during the solutiuon of the unsteady
        //  adjoint equations
        incompressibleVars(const incompressibleVars& vs);

        //- Clone the incompressibleVars
        virtual autoPtr<variablesSet> clone() const;



    //- Destructor
    virtual ~incompressibleVars() = default;


    // Member Functions

        // Access

            // Access to fields that will be later used for the solution of
            // the adjoint equations. Might be averaged or not, depending on
            // the corresponding switch
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                //- Return const reference to pressure
                const volScalarField& p() const;

                //- Return reference to pressure
                volScalarField& p();

                //- Return const reference to velocity
                const volVectorField& U() const;

                //- Return reference to velocity
                volVectorField& U();

                //- Return const reference to volume flux
                const surfaceScalarField& phi() const;

                //- Return reference to volume flux
                surfaceScalarField& phi();


            // Access to instantaneous fields. Solvers should use these fields
            // to execute a solver iteration
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                //- Return const reference to pressure
                inline const volScalarField& pInst() const;

                //- Return reference to pressure
                inline volScalarField& pInst();

                //- Return const reference to velocity
                inline const volVectorField& UInst() const;

                //- Return reference to velocity
                inline volVectorField& UInst();

                //- Return const reference to volume flux
                inline const surfaceScalarField& phiInst() const;

                //- Return reference to volume flux
                inline surfaceScalarField& phiInst();


            //- Return const reference to transport model
            inline const singlePhaseTransportModel& laminarTransport() const;

            //- Return reference to transport model
            inline singlePhaseTransportModel& laminarTransport();

            //- Return const reference to the turbulence model
            inline const autoPtr<incompressible::turbulenceModel>&
                turbulence() const;

            //- Return  reference to the turbulence model
            inline autoPtr<incompressible::turbulenceModel>& turbulence();

            //- Return const reference to the turbulence model variables
            inline const autoPtr<incompressible::RASModelVariables>&
                RASModelVariables() const;

            //- Return reference to the turbulence model variables
            inline autoPtr<incompressible::RASModelVariables>&
                RASModelVariables();

            //- Restore field values to the initial ones
            void restoreInitValues();

            //- Reset mean fields to zero
            void resetMeanFields();

            //- Compute mean fields on the fly
            void computeMeanFields();

            //- Compute mean fields on the fly: unsteady runs
            void computeMeanUnsteadyFields();

            //- correct boundaryconditions for all volFields
            void correctBoundaryConditions();

            //- Return storeInitValues bool
            bool storeInitValues() const;

            //- Return computeMeanFields bool
            bool computeMeanFields() const;

            /*
            //- Return const reference to the temperature
            const volScalarField& T() const;

            //- Return reference to the temperature
            volScalarField& T();
            */

            //- Transfer the fields of another variablesSet to this
            virtual void transfer(variablesSet& vars);

            //- Update nut
            virtual void validateTurbulence();

            //- Write dummy turbulent fields to allow for continuation in
            //- multi-point, turbulent runs
            bool write() const;

            //- Return const reference to the fvMesh
            fvMesh& mesh() const;

            //- Return const reference to the solverDict
            const solverControl& solverControlReference() const;

            //- const reference to the wordList containing the names of
            //- the primal fields that are allocated.
            const wordList& allocatedFieldNames() const override;

            //- const reference to the solution directions
            //- of the solver.
            const Vector<label>& solDirs() const;

            //- Set averaging start time in solverControl.
            //  Required by unsteady runs
            void adjustAverageStartTime(const scalar& offset) override;

            //- Set writeOption of all primal fields
            void setWriteOption(IOobject::writeOption w);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "incompressibleVarsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
