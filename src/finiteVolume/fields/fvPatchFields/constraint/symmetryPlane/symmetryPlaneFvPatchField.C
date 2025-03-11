/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "symmetryPlaneFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(p, iF),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p))
{}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchField<Type>(ptf, p, iF, mapper),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p))
{
    if (!isType<symmetryPlaneFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchField<Type>(p, iF, dict),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p, dict))
{
    if (!isType<symmetryPlaneFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf
)
:
    basicSymmetryFvPatchField<Type>(ptf),
    symmetryPlanePatch_(ptf.symmetryPlanePatch_)
{}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(ptf, iF),
    symmetryPlanePatch_(ptf.symmetryPlanePatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::symmetryPlaneFvPatchField<Type>::snGrad() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type: treat like zero-gradient
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        const symmTensor rot(I - 2.0*sqr(symmetryPlanePatch_.n()));

        const Field<Type> pif(this->patchInternalField());

        const auto& dc = this->patch().deltaCoeffs();

        return
        (
            (0.5*dc)
          * (transform(rot, pif) - pif)
        );
    }
}


template<class Type>
void Foam::symmetryPlaneFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type: treat like zero-gradient
        this->extrapolateInternal();
    }
    else
    {
        const symmTensor rot(I - 2.0*sqr(symmetryPlanePatch_.n()));

        const Field<Type> pif(this->patchInternalField());

        Field<Type>::operator=(0.5*(pif + transform(rot, pif)));
    }

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::symmetryPlaneFvPatchField<Type>::snGradTransformDiag() const
{
    const vector diag(cmptMag(symmetryPlanePatch_.n()));

    return tmp<Field<Type>>::New
    (
        this->size(),
        transformMask<Type>
        (
            //pow<vector, pTraits<Type>::rank>(diag)
            pow
            (
                diag,
                pTraits
                <
                    typename powProduct<vector, pTraits<Type>::rank>::type
                >::zero
            )
        )
    );
}


// ************************************************************************* //
