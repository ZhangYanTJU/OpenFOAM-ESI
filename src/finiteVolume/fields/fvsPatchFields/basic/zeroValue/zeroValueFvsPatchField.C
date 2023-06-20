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

#include "zeroValueFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::zeroValueFvsPatchField<Type>::zeroValueFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFvsPatchField<Type>::zeroValueFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{
    fvsPatchFieldBase::readDict(dict);
}


template<class Type>
Foam::zeroValueFvsPatchField<Type>::zeroValueFvsPatchField
(
    const zeroValueFvsPatchField<Type>& pfld,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper&
)
:
    // Field is zero. No mapping
    parent_bctype(pfld, p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFvsPatchField<Type>::zeroValueFvsPatchField
(
    const zeroValueFvsPatchField<Type>& pfld,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    // Field is zero
    parent_bctype(pfld, pfld.patch(), iF, Type(Zero))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvsPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // No contribution from internal values
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvsPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Patch field is zero
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvsPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroValueFvsPatchField<Type>::gradientBoundaryCoeffs() const
{
    // Patch field is zero
    return tmp<Field<Type>>::New(this->size(), Foam::zero{});
}


template<class Type>
void Foam::zeroValueFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    // Without writeValueEntry() since the value == zero
}


// ************************************************************************* //
