/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "basicSymmetryFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict)
{
    this->evaluate();
}


template<class Type>
Foam::basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFvPatchField<Type>::snGrad() const
{
    if (pTraits<Type>::rank == 0)
    {
        // Transform-invariant types
        return tmp<Field<Type>>::New(this->size(), Zero);
    }
    else
    {
        symmTensorField rot(I - 2.0*sqr(this->patch().nf()));

        Field<Type> pif(this->patchInternalField());

        return
        (
            (transform(rot, pif) - pif)
          * (this->patch().deltaCoeffs()/2.0)
        );
    }
}


template<class Type>
void Foam::basicSymmetryFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if (pTraits<Type>::rank == 0)
    {
        // Transform-invariant types
        Field<Type>::operator=(this->patchInternalField());
    }
    else
    {
        symmTensorField rot(I - 2.0*sqr(this->patch().nf()));

        Field<Type> pif(this->patchInternalField());

        Field<Type>::operator=((pif + transform(rot, pif))/2.0);
    }

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFvPatchField<Type>::snGradTransformDiag() const
{
    tmp<vectorField> diag(cmptMag(this->patch().nf()));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// ************************************************************************* //
