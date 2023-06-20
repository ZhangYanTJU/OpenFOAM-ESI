/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

#include "zeroValueFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::zeroValueFvPatchField<Type>::zeroValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFvPatchField<Type>::zeroValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{
    fvPatchFieldBase::readDict(dict);
}


template<class Type>
Foam::zeroValueFvPatchField<Type>::zeroValueFvPatchField
(
    const zeroValueFvPatchField<Type>& pfld,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    // Field is zero. No mapping
    parent_bctype(pfld, p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFvPatchField<Type>::zeroValueFvPatchField
(
    const zeroValueFvPatchField<Type>& pfld,
    const DimensionedField<Type, volMesh>& iF
)
:
    // Field is zero
    parent_bctype(pfld, pfld.patch(), iF, Type(Zero))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // No contribution from internal values
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Patch field is zero
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    // Patch field is zero
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
void Foam::zeroValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    // Without writeValueEntry() since the value == zero
}


// ************************************************************************* //
