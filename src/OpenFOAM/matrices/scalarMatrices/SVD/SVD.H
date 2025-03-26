/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::SVD

Description
    Singular value decomposition of a rectangular matrix.

SourceFiles
    SVD.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SVD_H
#define Foam_SVD_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Cmpt> class Tensor;

/*---------------------------------------------------------------------------*\
                      Class SVD Declaration
\*---------------------------------------------------------------------------*/

class SVD
{
    // Private Data

        //- Rectangular matrix with the same dimensions as the input
        scalarRectangularMatrix U_;

        //- Square matrix V
        scalarRectangularMatrix V_;

        //- The singular values
        DiagonalMatrix<scalar> S_;

        //- Convergence flag
        bool converged_;

        //- The number of zero singular values
        label nZeros_;


public:

    // Generated Methods

        //- No copy construct
        SVD(const SVD&) = delete;

        //- No copy assignment
        void operator=(const SVD&) = delete;


    // Constructors

        //- Construct from a rectangular Matrix
        explicit SVD
        (
            const scalarRectangularMatrix& A,
            const scalar minCondition = 0
        );


    // Access Functions

        //- Return U
        const scalarRectangularMatrix& U() const noexcept { return U_; }

        //- Return the square matrix V
        const scalarRectangularMatrix& V() const noexcept { return V_; }

        //- Return the singular values
        const scalarDiagonalMatrix& S() const noexcept { return S_; }

        //- Return the minimum non-zero singular value
        bool converged() const noexcept { return converged_; }

        //- Return the number of zero singular values
        label nZeros() const noexcept { return nZeros_; }

        //- Return the minimum non-zero singular value
        scalar minNonZeroS() const;


    // Member Functions

        //- Return the matrix product V S^(-1) U^T (the pseudo inverse)
        scalarRectangularMatrix VSinvUt() const;


    // Static Member Functions

        //- Return the pseudo inverse of the given matrix
        static scalarRectangularMatrix pinv
        (
            const scalarRectangularMatrix& A,
            const scalar minCondition = 0
        );

        //- Return the pseudo inverse of the given tensor
        static Tensor<scalar> pinv
        (
            const Tensor<scalar>& A,
            const scalar minCondition = 0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
