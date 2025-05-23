/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "uniformJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformJumpAMIFvPatchField<Type>::uniformJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<Type>(p, iF),
    jumpTable_()
{}


template<class Type>
Foam::uniformJumpAMIFvPatchField<Type>::uniformJumpAMIFvPatchField
(
    const uniformJumpAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<Type>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_.clone())
{}


template<class Type>
Foam::uniformJumpAMIFvPatchField<Type>::uniformJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<Type>(p, iF),
    jumpTable_()
{
    if (this->cyclicAMIPatch().owner())
    {
        jumpTable_ = Function1<Type>::New("jumpTable", dict, &this->db());
    }

    if (!this->readValueEntry(dict))
    {
        this->evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::uniformJumpAMIFvPatchField<Type>::uniformJumpAMIFvPatchField
(
    const uniformJumpAMIFvPatchField<Type>& ptf
)
:
    fixedJumpAMIFvPatchField<Type>(ptf),
    jumpTable_(ptf.jumpTable_.clone())
{}


template<class Type>
Foam::uniformJumpAMIFvPatchField<Type>::uniformJumpAMIFvPatchField
(
    const uniformJumpAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<Type>(ptf, iF),
    jumpTable_(ptf.jumpTable_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformJumpAMIFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicAMIPatch().owner())
    {
        this->jump_ = jumpTable_->value(this->db().time().value());
    }

    fixedJumpAMIFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformJumpAMIFvPatchField<Type>::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<Type>::write(os);
    if (this->cyclicAMIPatch().owner())
    {
        jumpTable_->writeData(os);
    }
}


// ************************************************************************* //
