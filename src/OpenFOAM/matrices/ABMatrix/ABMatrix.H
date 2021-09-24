/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::ABMatrix

Description

Usage

Notes

See also
    Test-ABMatrix.C

SourceFiles
    ABMatrix.C
    ABMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef ABMatrix_H
#define ABMatrix_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class ABMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class ABMatrix
{
public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> SMatrix;
    typedef RectangularMatrix<cmptType> RMatrix;


    //- Options for the output matrix forms of QRMatrix
    enum modes : uint8_t
    {
        FULL = 1,         //!<
        ECONOMY = 2,      //!<
    };

    //-
    /* enum inplace : bool
    {
        FALSE = false,      //!<
        TRUE  = true        //!<
    }; */

    //- Options for the computation of column pivoting
    enum pivoting : bool
    {
        FALSE = false,      //!< switches off column pivoting
        TRUE  = true        //!< switches on column pivoting
    };


private:

    // Private Data

        //-
        const modes mode_;

        //-
        const pivoting pivot_;

        //-
        label k_;

        //-
        MatrixType QR_;

        //-
        List<cmptType> Rdiag_;

        //- Permutation list (if column-pivoting is on)
        List<label> P_;


    // Private Member Functions

        //-
        void mode(const MatrixType& A);

        //-
        void pivot(MatrixType& A);


public:

    // Constructors

        //- Construct null
        ABMatrix();

        //- Construct QRMatrix without performing the decomposition
        ABMatrix
        (
            const modes mode,
            const pivoting pivot = pivoting::FALSE
        );

        //- Construct ABMatrix and perform the QR decomposition
        ABMatrix
        (
            MatrixType& A,
            const modes mode = modes::ECONOMY,
            const pivoting pivot = pivoting::FALSE
        );

        //- Construct ABMatrix and perform the QR decomposition
        ABMatrix
        (
            const MatrixType& A,
            const modes mode = modes::ECONOMY,
            const pivoting pivot = pivoting::FALSE
        );


    // Member Functions

        // Access

            //-
            inline const MatrixType& QR() const noexcept;

            //-
            inline const List<cmptType>& Rdiag() const noexcept;

            //-
            inline const List<label>& P() const noexcept;

            //- Create and return the permutation matrix
            inline SMatrix permut() const;


        // Evaluate

            //- Return the unitary similarity matrix
            //  Includes implicit round-to-zero as mutable operation
            MatrixType Q() const;

            //-
            MatrixType Q(const MatrixType& A) const;

            //- Return the upper triangular matrix
            MatrixType R() const;

            //-
            MatrixType R(const MatrixType& A) const;

            //-
            MatrixType H() const;

            //-
            MatrixType H(const MatrixType& A) const;

            //- Compute QR decomposition according to constructor settings
            void decompose(MatrixType& A);

            //- Compute QR decomposition according to constructor settings
            void decompose(const MatrixType& A);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ABMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ABMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
