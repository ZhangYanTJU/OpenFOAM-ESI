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
    const basicSymmetryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::basicSymmetryFvPatchField<Type>::snGrad(UList<Type>& result) const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : treat like zero-gradient
        result = Foam::zero{};
    }
    else
    {
        // Get patch internal field, stored temporarily in result
        this->patchInternalField(result);
        const auto& pif = result;

        tmp<vectorField> tnHat = this->patch().nf();
        const auto& nHat = tnHat();

        const auto& dc = this->patch().deltaCoeffs();

        const label len = result.size();

        // (dc/2.0)*(transform(I - 2.0*sqr(nHat), iF) - iF);

        for (label i = 0; i < len; ++i)
        {
            result[i] =
            (
                (0.5*dc[i])
              * (transform(I - 2.0*sqr(nHat[i]), pif[i]) - pif[i])
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFvPatchField<Type>::snGrad() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : treat like zero-gradient
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        auto tresult = tmp<Field<Type>>::New(this->size());
        this->snGrad(static_cast<UList<Type>&>(tresult.ref()));
        return tresult;
    }
}


template<class Type>
void Foam::basicSymmetryFvPatchField<Type>::evaluate(const Pstream::commsTypes)
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
        tmp<vectorField> nHat = this->patch().nf();

        const Field<Type> pif(this->patchInternalField());

        Field<Type>::operator=
        (
            0.5*(pif + transform(I - 2.0*sqr(nHat), pif))
        );
    }

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::basicSymmetryFvPatchField<Type>::snGradTransformDiag() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type
        // FatalErrorInFunction
        //     << "Should not be called for this type"
        //     << ::Foam::abort(FatalError);
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        tmp<vectorField> diag(cmptMag(this->patch().nf()));

        return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
    }
}


// ************************************************************************* //
