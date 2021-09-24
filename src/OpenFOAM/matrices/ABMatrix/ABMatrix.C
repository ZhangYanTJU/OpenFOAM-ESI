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

\*---------------------------------------------------------------------------*/

#include "ABMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
void Foam::ABMatrix<MatrixType>::mode(const MatrixType& A)
{
    if ((mode_ == modes::ECONOMY) && (A.m() <= A.n()))
    {
        FatalErrorInFunction
            << "Economy QR decomposition is by definition limited to matrices "
            << "wherein rows > columns. Please use FULL decomposition, instead."
            << abort(FatalError);
    }

    if (mode_ == modes::FULL)
    {
        k_ = A.m();
    }
    else if (mode_ == modes::ECONOMY)
    {
        k_ = min(A.m(), A.n());
    }
}


template<class MatrixType>
void Foam::ABMatrix<MatrixType>::pivot
(
    MatrixType& A
)
{
    const label nCols = A.n();
    const label nIter = min(A.m() - 1, A.n());

    // Initialise permutation vector, and column norm vector
    P_ = identity(nCols);

    // Initialise vector norms of each column of A
    List<scalar> colNorms(nCols);
    for (label k = 0; k < nCols; ++k)
    {
        colNorms[k] = A.columnNorm(k, true);
    }

    // Loop through all subcolumns of which the diagonal elem is the first elem
    for (label k = 0; k < nIter; ++k)
    {
        const labelRange colRange(k, nCols);
        const SubList<scalar> subColNorms(colNorms, colRange);

        // Column pivoting
        const label maxColNormi =
            std::distance
            (
                subColNorms.cbegin(),
                std::max_element(subColNorms.cbegin(), subColNorms.cend())
            );

        // Swap R_, P_ and colNorms_ according to pivot column if the current
        // column is not the max norm column by avoiding any column swap where
        // the leading elem is very small
        if (maxColNormi != 0 && SMALL < mag(A(k, k + maxColNormi)))
        {
            const RMatrix R1(A.subColumn(k));
            const RMatrix R2(A.subColumn(maxColNormi + k));
            A.subColumn(k) = R2;
            A.subColumn(maxColNormi + k) = R1;

            Swap(P_[k], P_[maxColNormi + k]);
            Swap(colNorms[k], colNorms[maxColNormi + k]);
        }

        // Update norms
        if (k < nIter - 1)
        {
            label q = k + 1;
            for (const auto& val : RMatrix(A.subRow(k, q)))
            {
                colNorms[q] -= magSqr(val);
                ++q;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::ABMatrix<MatrixType>::ABMatrix()
:
    mode_(modes::ECONOMY),
    pivot_(pivoting::FALSE),
    k_(0),
    QR_(),
    Rdiag_(),
    P_()
{}


template<class MatrixType>
Foam::ABMatrix<MatrixType>::ABMatrix
(
    const modes mode,
    const pivoting pivot
)
:
    mode_(mode),
    pivot_(pivot),
    k_(0),
    QR_(),
    Rdiag_(),
    P_()
{}


template<class MatrixType>
Foam::ABMatrix<MatrixType>::ABMatrix
(
    MatrixType& A,
    const modes mode,
    const pivoting pivot
)
:
    mode_(mode),
    pivot_(pivot),
    k_(0),
    QR_(),
    Rdiag_(),
    P_()
{
    decompose(A);
}


template<class MatrixType>
Foam::ABMatrix<MatrixType>::ABMatrix
(
    const MatrixType& A,
    const modes mode,
    const pivoting pivot
)
:
    mode_(mode),
    pivot_(pivot),
    k_(0),
    QR_(),
    Rdiag_(),
    P_()
{
    decompose(A);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MatrixType>
void Foam::ABMatrix<MatrixType>::decompose
(
    MatrixType& A
)
{
    Info<< "inline" << endl;

    mode(A);

    if (pivot_)
    {
        Info<< "pivoting: on" << endl;
        pivot(A);
    }

    const label m = A.m();
    const label n = A.n();

    Rdiag_ = List<cmptType>(n, Zero);

    for (label k = 0; k < n; ++k)
    {
        // Compute 2-norm of k-th column without under/overflow
        scalar nrm = scalar(0);

        for (label i = k; i < m; ++i)
        {
            nrm = hypot(nrm, A[i][k]);
        }

        if (nrm != scalar(0))
        {
            // Form k-th Householder vector
            if (A[k][k] < 0)
            {
                nrm = -nrm;
            }

            for (label i = k; i < m; ++i)
            {
                A[i][k] /= nrm;
            }

            A[k][k] += scalar(1);

            // Apply transformation to remaining columns
            for (label j = k + 1; j < n; ++j)
            {
                scalar s = 0;

                for (label i = k; i < m; ++i)
                {
                    s += A[i][k]*A[i][j];
                }

                s /= -A[k][k];

                for (label i = k; i < m; ++i)
                {
                    A[i][j] += s*A[i][k];
                }
            }
        }

        Rdiag_[k] = -nrm;
    }
}


template<class MatrixType>
void Foam::ABMatrix<MatrixType>::decompose
(
    const MatrixType& A
)
{
    Info<< "offline" << endl;
    QR_ = A;

    decompose(QR_);
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::Q
(
    const MatrixType& A
) const
{
    const label m = A.m();

    MatrixType Q(m, k_, Zero);

    for (label k = k_ - 1; k >= 0; --k)
    {
        for (label i = 0; i < m; ++i)
        {
            Q[i][k] = scalar(0);
        }

        Q[k][k] = scalar(1);

        for (label j = k; j < k_; ++j)
        {
            if (A[k][k] != scalar(0))
            {
                scalar s = 0;

                for (label i = k; i < m; ++i)
                {
                    s += A[i][k]*Q[i][j];
                }

                s /= -A[k][k];

                for (label i = k; i < m; ++i)
                {
                    Q[i][j] += s*A[i][k];
                }
            }
        }
    }

    return Q;
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::Q() const
{
    return Q(QR_);
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::R
(
    const MatrixType& A
) const
{
    const label n = A.n();

    MatrixType R(k_, n, Zero);

    for (label i = 0; i < k_; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (i < j)
            {
                R[i][j] = A[i][j];
            }
            else if (i == j)
            {
                R[i][j] = Rdiag_[i];
            }
            else
            {
                R[i][j] = scalar(0);
            }
        }
    }

    return R;
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::R() const
{
    return R(QR_);
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::H
(
    const MatrixType& A
) const
{
    const label m = A.m();
    const label n = A.n();

    MatrixType H(m, n, Zero);

    for (label i = 0; i < m; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (i >= j)
            {
                H[i][j] = A[i][j];
            }
            else
            {
                H[i][j] = scalar(0);
            }
        }
    }

    return H;
}


template<class MatrixType>
MatrixType Foam::ABMatrix<MatrixType>::H() const
{
    return H(QR_);
}


// ************************************************************************* //
