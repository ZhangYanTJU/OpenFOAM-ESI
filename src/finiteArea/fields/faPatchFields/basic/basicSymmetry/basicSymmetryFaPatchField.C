/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

\*---------------------------------------------------------------------------*/

#include "basicSymmetryFaPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    this->evaluate();
}


template<class Type>
Foam::basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFaPatchField<Type>::snGrad() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : treat like zero-gradient
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        tmp<vectorField> nHat = this->patch().edgeNormals();

        const auto& dc = this->patch().deltaCoeffs();

        const Field<Type> pif(this->patchInternalField());

        return
        (
            (0.5*dc)
          * (transform(I - 2.0*sqr(nHat), pif) - pif)
        );
    }
}


template<class Type>
void Foam::basicSymmetryFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : treat like zero-gradient
        this->extrapolateInternal();
    }
    else
    {
        tmp<vectorField> nHat = this->patch().edgeNormals();

        const Field<Type> pif(this->patchInternalField());

        Field<Type>::operator=
        (
            0.5*(pif + transform(I - 2.0*sqr(nHat), pif))
        );
    }

    transformFaPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFaPatchField<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().edgeNormals()));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// ************************************************************************* //
