/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024-2025 OpenCFD Ltd.
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

#include "wedgeFvPatch.H"
#include "wedgeFvPatchField.H"
#include "transformField.H"
#include "symmTransform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF)
{}


template<class Type>
Foam::wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper)
{
    if (!isType<wedgeFvPatch>(this->patch()))
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
Foam::wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)  // "value" is NO_READ
{
    if (!isType<wedgeFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    evaluate();
}


template<class Type>
Foam::wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(ptf, iF)
{}


template<class Type>
Foam::wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf
)
:
    wedgeFvPatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::wedgeFvPatchField<Type>::snGrad() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : treat like zero-gradient
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        const auto& rot = refCast<const wedgeFvPatch>(this->patch()).cellT();

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
void Foam::wedgeFvPatchField<Type>::evaluate(const Pstream::commsTypes)
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
        const auto& rot = refCast<const wedgeFvPatch>(this->patch()).faceT();

        fvPatchField<Type>::operator==
        (
            transform(rot, this->patchInternalField())
        );
     }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::wedgeFvPatchField<Type>::snGradTransformDiag() const
{
    const auto& rot = refCast<const wedgeFvPatch>(this->patch()).cellT();

    const vector diagV = 0.5*(I - rot).diag();

    return tmp<Field<Type>>::New
    (
        this->size(),
        transformMask<Type>
        (
            pow
            (
                diagV,
                pTraits
                <
                    typename powProduct<vector, pTraits<Type>::rank>::type
                >::zero
            )
        )
    );
}


// ************************************************************************* //
