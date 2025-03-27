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
    Foam::simpleMatrix

Description
    A simple square matrix solver with scalar coefficients.

SourceFiles
    simpleMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_simpleMatrix_H
#define Foam_simpleMatrix_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class simpleMatrix;
template<class Type> Ostream& operator<<(Ostream&, const simpleMatrix<Type>&);

/*---------------------------------------------------------------------------*\
                           Class simpleMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class simpleMatrix
:
    public scalarSquareMatrix
{
    // Private Data

        Field<Type> source_;

public:

    // Constructors

        //- Construct given size
        //  Note: this does not initialise the coefficients or the source.
        explicit simpleMatrix(const label n);

        //- Construct given size and initial values for coefficients and source
        simpleMatrix(const label n, const scalar, const Type&);

        //- Construct given size, zero values for coefficients and source
        simpleMatrix(const label n, Foam::zero);

        //- Construct from components
        simpleMatrix(const scalarSquareMatrix&, const Field<Type>&);

        //- Copy construct
        simpleMatrix(const simpleMatrix<Type>&) = default;

        //- Construct from Istream
        explicit simpleMatrix(Istream&);


    // Member Functions

        //- Return access to the source
        Field<Type>& source() noexcept { return source_; }

        //- Return const-access to the source
        const Field<Type>& source() const noexcept { return source_; }


        //- Solve the matrix using Gaussian elimination with pivoting
        //- and return the solution
        Field<Type> solve() const;

        //- Solve the matrix using LU decomposition with pivoting
        //- and return the solution
        Field<Type> LUsolve() const;


    // Member Operators

        //- Copy assignment
        void operator=(const simpleMatrix<Type>&);


    // Ostream Operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const simpleMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

template<class Type>
simpleMatrix<Type> operator+
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator-
(
    const simpleMatrix<Type>&,
    const simpleMatrix<Type>&
);

template<class Type>
simpleMatrix<Type> operator*
(
    const scalar,
    const simpleMatrix<Type>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "simpleMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
