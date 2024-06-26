/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "emptyFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::emptyFvPatchField<Type>::emptyFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{}


template<class Type>
Foam::emptyFvPatchField<Type>::emptyFvPatchField
(
    const emptyFvPatchField<Type>&,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper&
)
:
    fvPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    if (!isType<emptyFvPatch>(p))
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
Foam::emptyFvPatchField<Type>::emptyFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    // Empty means empty, so no patchType override
    // with fvPatchFieldBase::readDict(dict);

    if (!isType<emptyFvPatch>(p))
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
Foam::emptyFvPatchField<Type>::emptyFvPatchField
(
    const emptyFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf.patch(), iF, Field<Type>())  // zero-sized
{}


template<class Type>
Foam::emptyFvPatchField<Type>::emptyFvPatchField
(
    const emptyFvPatchField<Type>& ptf
)
:
    emptyFvPatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::emptyFvPatchField<Type>::updateCoeffs()
{
    // Check moved to checkMesh.
    // Test here fails if multiple empty patches.
}


template<class Type>
void Foam::emptyFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    // Never write "value"
}


// ************************************************************************* //
