/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::LUscalarMatrix

Description
    Class to perform the LU decomposition on a symmetric matrix.

SourceFiles
    LUscalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LUscalarMatrix_H
#define Foam_LUscalarMatrix_H

#include "scalarMatrices.H"
#include "labelList.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduMatrix;
class procLduMatrix;

/*---------------------------------------------------------------------------*\
                           Class LUscalarMatrix Declaration
\*---------------------------------------------------------------------------*/

class LUscalarMatrix
:
    public scalarSquareMatrix
{
    // Private Data

        //- Communicator to use
        const label comm_;

        //- Processor matrix offsets
        labelList procOffsets_;

        //- The pivot indices used in the LU decomposition
        labelList pivotIndices_;


    // Private Member Functions

        //- Convert the given lduMatrix into this LUscalarMatrix
        void convert
        (
            const lduMatrix& ldum,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Convert the given list of procLduMatrix into this LUscalarMatrix
        //- on the master processor
        void convert(const PtrList<procLduMatrix>& lduMatrices);


        //- Print the ratio of the mag-sum of the off-diagonal coefficients
        //- to the mag-diagonal
        void printDiagonalDominance() const;


public:

    // Declare name of the class and its debug switch
    ClassName("LUscalarMatrix");


    // Constructors

        //- Default construct
        LUscalarMatrix() noexcept;

        //- Construct from and perform LU decomposition of the given matrix
        explicit LUscalarMatrix(const scalarSquareMatrix& mat);

        //- Construct from lduMatrix and perform LU decomposition.
        //- In parallel it assembles the matrix on the master.
        LUscalarMatrix
        (
            const lduMatrix& ldum,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Member Functions

        //- Perform the LU decomposition of the matrix
        void decompose(const scalarSquareMatrix& mat);

        //- Solve the linear system with the given source
        //  and returning the solution in the Field argument x.
        //  This function may be called with the same field for x and source.
        template<class Type>
        void solve(List<Type>& x, const UList<Type>& source) const;

        //- Solve the linear system with the given source
        //  returning the solution
        template<class Type>
        tmp<Field<Type>> solve(const UList<Type>& source) const;

        //- Set M to the inverse of this square matrix
        void inv(scalarSquareMatrix& M) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LUscalarMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
