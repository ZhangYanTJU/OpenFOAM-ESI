/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

Description
    Various tools for test applications.

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "complex.H"
#include "Matrix.H"
#include "Random.H"
#include <chrono>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total number of unit tests
unsigned nTest_ = 0;


// Total number of failed unit tests
unsigned nFail_ = 0;


// Return the absolute tolerance value for bitwise comparisons of floatScalars
floatScalar getTol(floatScalar)
{
    return 1e-2;
}


// Return the absolute tolerance value for bitwise comparisons of doubleScalars
doubleScalar getTol(doubleScalar)
{
    return 1e-10;
}


// Return the absolute tolerance value for bitwise comparisons of doubleScalars
doubleScalar getTol(complex)
{
    return 1e-10;
}


// Create a random Matrix (real or complex)
template<class MatrixType>
MatrixType makeRandomMatrix
(
    const labelPair& dims,
    Random& rndGen
)
{
    MatrixType mat(dims);

    if constexpr (std::is_same_v<complex, typename MatrixType::cmptType>)
    {
        for (auto& x : mat)
        {
            x.real(rndGen.GaussNormal<scalar>());
            x.imag(rndGen.GaussNormal<scalar>());
        }
    }
    else
    {
        std::generate
        (
            mat.begin(),
            mat.end(),
            [&]{ return rndGen.GaussNormal<scalar>(); }
        );
    }

    return mat;
}


// Copy an initializer list into a DiagonalMatrix
template<class Type>
void assignMatrix
(
    UList<Type>& A,
    std::initializer_list<Type> list
)
{
    std::copy(list.begin(), list.end(), A.begin());
}


// Copy an initializer list into a SymmetricSquareMatrix.
template<class Form, class Type>
void assignMatrix
(
    Matrix<Form, Type>& A,
    std::initializer_list<typename Matrix<Form, Type>::cmptType> list
)
{
    const label nargs = list.size();

    if (nargs != A.size())
    {
        FatalErrorInFunction
            << "Mismatch in matrix dimension ("
            << A.m() << ", "
            << A.n() << ") and number of args (" << nargs << ')' << nl
            << exit(FatalError);
     }

    std::copy(list.begin(), list.end(), A.begin());
}


// Return a copy of the Matrix collapsed into one dimension.
template<class Form, class Type>
List<Type> flt
(
    const Matrix<Form, Type>& A,
    const bool rowMajorOrder = true
)
{
   List<Type> flatMatrix(A.size());

   if (rowMajorOrder)
   {
       std::copy(A.cbegin(), A.cend(), flatMatrix.begin());
   }
   else
   {
       for (label j = 0; j < A.n(); ++j)
       {
           for (label i = 0; i < A.m(); ++i)
           {
               flatMatrix[i + j*A.m()] = A(i, j);
           }
       }
   }

   return flatMatrix;
}


// Compare two values or two containers (elementwise), and print output.
// Do ++nFail_ if values of two objects are not equal within a given tolerance.
// The function is converted from PEP-485.
template<class Type1, class Type2 = Type1>
void cmp
(
    const word& msg,
    const Type1& x,
    const Type2& y,
    const scalar absTol = 0,    //<! useful for cmps near zero
    const scalar relTol = 1e-8, //<! are values the same within 8 decimals
    const bool verbose = false
)
{
    const auto notEqual = [=](const auto& a, const auto& b) -> bool
    {
        return
        (
            Foam::max(absTol, relTol*Foam::max(Foam::mag(a), Foam::mag(b)))
          < Foam::mag(a - b)
        );
    };

    unsigned nFail = 0;

    if constexpr
    (
        std::is_floating_point_v<Type1> || std::is_same_v<complex, Type1>
    )
    {
        if (notEqual(x, y))
        {
            ++nFail;
        }
    }
    else
    {
        for (label i = 0; i < x.size(); ++i)
        {
            if (notEqual(x[i], y[i]))
            {
                ++nFail;
            }
        }
    }

    if (verbose)
    {
        Info<< msg << x << "?=" << y << endl;
    }

    if (nFail)
    {
        Info<< nl
            << "        #### Fail in " << nFail << " comps ####" << nl << endl;
        ++nFail_;
    }
    ++nTest_;
}


// Compare two Booleans, and print output.
// Do ++nFail_ if two Booleans are not equal.
void cmp
(
    const word& msg,
    const bool x,
    const bool y,
    const bool verbose = false
)
{
    unsigned nFail = 0;

    if (x != y)
    {
        ++nFail;
    }

    if (verbose)
    {
        Info<< msg << x << "?=" << y << endl;
    }

    if (nFail)
    {
        Info<< nl
            << "        #### Fail in " << nFail << " comps ####" << nl << endl;
        ++nFail_;
    }
    ++nTest_;
}


// Print OpenFOAM matrix in NumPy format
template<class MatrixType>
void InfoNumPyFormat
(
    const MatrixType& mat,
    const word objName
)
{
    Info<< objName << ": " << mat.m() << "x" << mat.n() << nl;
    for (label m = 0; m < mat.m(); ++m)
    {
        Info<< "[";
        for (label n = 0; n < mat.n(); ++n)
        {
            if (n == mat.n() - 1)
            {
                Info<< mat(m,n);
            }
            else
            {
                Info<< mat(m,n) << ",";
            }
        }
        if (m == mat.m() - 1)
        {
            Info << "]" << nl;
        }
        else
        {
            Info << "]," << nl;
        }
    }
}


// ************************************************************************* //
