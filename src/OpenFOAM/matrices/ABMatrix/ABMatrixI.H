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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MatrixType>
inline const MatrixType&
Foam::ABMatrix<MatrixType>::QR() const noexcept
{
    return QR_;
}


template<class MatrixType>
inline const Foam::List<typename Foam::ABMatrix<MatrixType>::cmptType>&
Foam::ABMatrix<MatrixType>::Rdiag() const noexcept
{
    return Rdiag_;
}


template<class MatrixType>
inline const Foam::List<Foam::label>&
Foam::ABMatrix<MatrixType>::P() const noexcept
{
    return P_;
}


template<class MatrixType>
inline Foam::SquareMatrix<typename Foam::ABMatrix<MatrixType>::cmptType>
Foam::ABMatrix<MatrixType>::permut() const
{
    SquareMatrix<cmptType> permut(P_.size(), Zero);

    forAll(P_, jcol)
    {
        permut(P_[jcol], jcol) = pTraits<cmptType>::one;
    }

    return permut;
}


// ************************************************************************* //
