/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::DICSmoother

Group
    grpLduMatrixSmoothers

Description
    Simplified diagonal-based incomplete Cholesky smoother for symmetric
    matrices.

    To improve efficiency, the residual is evaluated after every nSweeps
    sweeps.

SourceFiles
    DICSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef DICSmoother_H
#define DICSmoother_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DICSmoother Declaration
\*---------------------------------------------------------------------------*/

class DICSmoother
:
    public lduMatrix::smoother
{
    // Private data

        //- The reciprocal preconditioned diagonal
        solveScalarField rD_;


public:

    //- Runtime type information
    TypeName("DIC");


    // Constructors

        //- Construct from matrix components
        DICSmoother
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    // Member Functions

        //- Smooth the solution for a given number of sweeps
        void smooth
        (
            solveScalarField& psi,
            const scalarField& source,
            const direction cmpt,
            const label nSweeps
        ) const;

        //- Smooth the solution for a given number of sweeps
        void scalarSmooth
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
