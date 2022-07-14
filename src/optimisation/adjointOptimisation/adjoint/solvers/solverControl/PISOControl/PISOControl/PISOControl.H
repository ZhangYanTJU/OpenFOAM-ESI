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
    Foam::PISOControl

Description
    PISO control class to supply convergence information/checks for
    the PISO loop.

\*---------------------------------------------------------------------------*/

#ifndef PISOControl_H
#define PISOControl_H

#include "solverControl.H"
#include "pisoControl.H"
#include "unsteadyTimeManipulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PISOControl Declaration
\*---------------------------------------------------------------------------*/

class PISOControl
:
    public solverControl,
    public pisoControl
{

protected:

    // Protected Data

        //- Optimisation type
        const word& managerType_;

        //- Pressure reference cell
        label pRefCell_;
        //
        //- Pressure reference value
        scalar pRefValue_;

        //- Time manipulation during unsteady optimisation
        unsteadyTimeManipulation& timeManip_;


private:

        //- No copy construct
        PISOControl(const PISOControl&) = delete;

        //- No copy assignment
        void operator=(const PISOControl&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("PISOControl");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            PISOControl,
            dictionary,
            (
                fvMesh& mesh,
                const word& managerType,
                const solver& solver
            ),
            (mesh, managerType, solver)
        );


    // Constructors

        //- Construct from mesh
        PISOControl
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    //- Destructor
    virtual ~PISOControl() = default;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<PISOControl> New
        (
            fvMesh& mesh,
            const word& managerType,
            const solver& solver
        );


    // Member Functions

        virtual bool read();

        // Access

            //- Return the solution dictionary
            inline virtual const dictionary dict() const;

            //- Return pressure reference cell
            inline label& pRefCell();

            //- Return pressure reference value
            inline scalar& pRefValue();

            //- Return the time span of each optimisation cycle
            inline scalar timeSpan() const;

            //- Update the start time of averaging for the next opt. cycle
            void updateAverageStartTime();


        // Solution control

            //- Whether to call time.write() or not
            virtual void writeFields() = 0;

            //- Use averaged fields?
            //  For solving the adjoint equations or computing sensitivities
            //  based on averaged fields
            inline bool useAveragedFields() const;


        // Evolution

            //- Loop
            virtual bool loop() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "PISOControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
