/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2015 Hrvoje Jasak, Wikki Ltd.
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
    gmresSolver

Group
    grpLduMatrixSolvers

Description
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

SourceFiles
    gmresSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gmresSolver_H
#define Foam_gmresSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gmresSolver Declaration
\*---------------------------------------------------------------------------*/

class gmresSolver
:
    public lduMatrix::solver
{
    // Private Data

        //- Preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr_;

        //- Krylov space dimension
        label nDirs_;


        //- Givens rotation
        void givensRotation
        (
            const solveScalar& H,
            const solveScalar& beta,
            solveScalar& c,
            solveScalar& s
        ) const;


public:

    //- Runtime type information
    TypeName("GMRES");


    // Generated Methods

        //- No copy construct
        gmresSolver(const gmresSolver&) = delete;

        //- No copy assignment
        void operator=(const gmresSolver&) = delete;


    // Constructors

        //- Construct from matrix components and solver data stream
        gmresSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~gmresSolver() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt=0
        ) const;

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
