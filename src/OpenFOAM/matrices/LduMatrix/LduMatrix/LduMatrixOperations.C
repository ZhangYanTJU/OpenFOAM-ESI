/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumDiag()
{
    const Field<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const Field<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();
    Field<DType>& Diag = diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face=0; face<l.size(); face++)
    {
        Diag[l[face]] += Lower[face];
        Diag[u[face]] += Upper[face];
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::negSumDiag()
{
    const Field<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const Field<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();
    Field<DType>& Diag = diag();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face=0; face<l.size(); face++)
    {
        Diag[l[face]] -= Lower[face];
        Diag[u[face]] -= Upper[face];
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::sumMagOffDiag
(
    Field<LUType>& sumOff
) const
{
    const Field<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const Field<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();

    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    for (label face = 0; face < l.size(); face++)
    {
        sumOff[u[face]] += cmptMag(Lower[face]);
        sumOff[l[face]] += cmptMag(Upper[face]);
    }
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::Field<Type>>
Foam::LduMatrix<Type, DType, LUType>::H(const Field<Type>& psi) const
{
    auto tHpsi = tmp<Field<Type>>::New(lduAddr().size(), Foam::zero{});

    if (hasLower() || hasUpper())
    {
        Type* __restrict__ HpsiPtr = tHpsi.ref().begin();

        const Type* __restrict__ psiPtr = psi.begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const LUType* __restrict__ lowerPtr = lower().begin();
        const LUType* __restrict__ upperPtr = upper().begin();

        const label nFaces = upper().size();

        for (label face=0; face<nFaces; face++)
        {
            HpsiPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
            HpsiPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
        }
    }

    return tHpsi;
}

template<class Type, class DType, class LUType>
Foam::tmp<Foam::Field<Type>>
Foam::LduMatrix<Type, DType, LUType>::H(const tmp<Field<Type>>& tpsi) const
{
    tmp<Field<Type>> tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::Field<Type>>
Foam::LduMatrix<Type, DType, LUType>::faceH(const Field<Type>& psi) const
{
    const Field<LUType>& Lower = const_cast<const LduMatrix&>(*this).lower();
    const Field<LUType>& Upper = const_cast<const LduMatrix&>(*this).upper();

    // Take references to addressing
    const labelUList& l = lduAddr().lowerAddr();
    const labelUList& u = lduAddr().upperAddr();

    auto tfaceHpsi = tmp<Field<Type>>::New(Lower.size());
    auto& faceHpsi = tfaceHpsi.ref();

    for (label face=0; face<l.size(); face++)
    {
        faceHpsi[face] = Upper[face]*psi[u[face]] - Lower[face]*psi[l[face]];
    }

    return tfaceHpsi;
}


template<class Type, class DType, class LUType>
Foam::tmp<Foam::Field<Type>>
Foam::LduMatrix<Type, DType, LUType>::faceH(const tmp<Field<Type>>& tpsi) const
{
    tmp<Field<Type>> tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator=(const LduMatrix& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    if (A.hasDiag())
    {
        diag() = A.diag();
    }

    if (A.hasUpper())
    {
        upper() = A.upper();
    }
    else
    {
        upperPtr_.reset(nullptr);
    }

    if (A.hasLower())
    {
        lower() = A.lower();
    }
    else
    {
        lowerPtr_.reset(nullptr);
    }

    if (A.hasSource())
    {
        source() = A.source();
    }

    interfacesUpper_ = A.interfacesUpper_;
    interfacesLower_ = A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator=(LduMatrix&& A)
{
    if (this == &A)
    {
        return;  // Self-assignment is a no-op
    }

    diagPtr_ = std::move(A.diagPtr_);
    upperPtr_ = std::move(A.upperPtr_);
    lowerPtr_ = std::move(A.lowerPtr_);
    sourcePtr_ = std::move(A.sourcePtr_);

    interfacesUpper_ = std::move(A.interfacesUpper_);
    interfacesLower_ = std::move(A.interfacesLower_);
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::negate()
{
    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (sourcePtr_)
    {
        sourcePtr_->negate();
    }

    negate(interfacesUpper_);
    negate(interfacesLower_);
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator+=(const LduMatrix& A)
{
    if (A.hasDiag())
    {
        diag() += A.diag();
    }

    if (A.hasSource())
    {
        source() += A.source();
    }

    if (symmetric() && A.symmetric())
    {
        upper() += A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {
            upper() = A.upper();
        }

        if (A.hasLower())
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorInFunction
            << "Unknown matrix type combination"
            << abort(FatalError);
    }

    interfacesUpper_ += A.interfacesUpper_;
    interfacesLower_ += A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator-=(const LduMatrix& A)
{
    if (A.hasDiag())
    {
        diag() -= A.diag();
    }

    if (A.hasSource())
    {
        source() -= A.source();
    }

    if (symmetric() && A.symmetric())
    {
        upper() -= A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.hasUpper())
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.hasUpper())
        {
            upper() = -A.upper();
        }

        if (A.hasLower())
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorInFunction
            << "Unknown matrix type combination"
            << abort(FatalError);
    }

    interfacesUpper_ -= A.interfacesUpper_;
    interfacesLower_ -= A.interfacesLower_;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator*=
(
    const scalarField& sf
)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (sourcePtr_)
    {
        *sourcePtr_ *= sf;
    }

    // Non-uniform scaling causes a symmetric matrix
    // to become asymmetric
    if (symmetric() || asymmetric())
    {
        Field<LUType>& upper = this->upper();
        Field<LUType>& lower = this->lower();

        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();

        for (label face=0; face<upper.size(); face++)
        {
            upper[face] *= sf[l[face]];
        }

        for (label face=0; face<lower.size(); face++)
        {
            lower[face] *= sf[u[face]];
        }
    }

    FatalErrorInFunction
        << "Scaling a matrix by scalarField is not currently supported\n"
           "because scaling interfacesUpper_ and interfacesLower_ "
           "require special transfers"
        << abort(FatalError);

    //interfacesUpper_ *= ;
    //interfacesLower_ *= sf;
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (sourcePtr_)
    {
        *sourcePtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }

    interfacesUpper_ *= s;
    interfacesLower_ *= s;
}


// ************************************************************************* //
