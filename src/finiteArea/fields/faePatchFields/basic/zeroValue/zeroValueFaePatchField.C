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

#include "zeroValueFaePatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::zeroValueFaePatchField<Type>::zeroValueFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFaePatchField<Type>::zeroValueFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    // Field is zero
    parent_bctype(p, iF, Type(Zero))
{
    faePatchFieldBase::readDict(dict);
}


template<class Type>
Foam::zeroValueFaePatchField<Type>::zeroValueFaePatchField
(
    const zeroValueFaePatchField<Type>& pfld,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper&
)
:
    // Field is zero. No mapping
    parent_bctype(pfld, p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroValueFaePatchField<Type>::zeroValueFaePatchField
(
    const zeroValueFaePatchField<Type>& pfld,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    // Field is zero
    parent_bctype(pfld, pfld.patch(), iF, Type(Zero))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::zeroValueFaePatchField<Type>::write(Ostream& os) const
{
    faePatchField<Type>::write(os);
    // Without writeValueEntry() since the value == zero
}


// ************************************************************************* //
