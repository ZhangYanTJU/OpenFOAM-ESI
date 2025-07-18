/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

#include "fixedMeanOutletInletFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedMeanOutletInletFvPatchField<Type>::fixedMeanOutletInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    outletInletFvPatchField<Type>(p, iF),
    meanValue_()
{}


template<class Type>
Foam::fixedMeanOutletInletFvPatchField<Type>::fixedMeanOutletInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    outletInletFvPatchField<Type>(p, iF),
    meanValue_(Function1<Type>::New("meanValue", dict, &this->db()))
{
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::fixedMeanOutletInletFvPatchField<Type>::fixedMeanOutletInletFvPatchField
(
    const fixedMeanOutletInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    outletInletFvPatchField<Type>(ptf, p, iF, mapper),
    meanValue_(ptf.meanValue_.clone())
{}


template<class Type>
Foam::fixedMeanOutletInletFvPatchField<Type>::fixedMeanOutletInletFvPatchField
(
    const fixedMeanOutletInletFvPatchField<Type>& ptf
)
:
    outletInletFvPatchField<Type>(ptf),
    meanValue_(ptf.meanValue_.clone())
{}


template<class Type>
Foam::fixedMeanOutletInletFvPatchField<Type>::fixedMeanOutletInletFvPatchField
(
    const fixedMeanOutletInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    outletInletFvPatchField<Type>(ptf, iF),
    meanValue_(ptf.meanValue_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedMeanOutletInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    Type meanValue = meanValue_->value(t);

    Field<Type> newValues(this->patchInternalField());

    Type meanValuePsi = gWeightedAverage(this->patch().magSf(), newValues);

    if (mag(meanValue) > SMALL && mag(meanValuePsi) > 0.5*mag(meanValue))
    {
        newValues *= mag(meanValue)/mag(meanValuePsi);
    }
    else
    {
        newValues += (meanValue - meanValuePsi);
    }

    this->refValue() = newValues;

    outletInletFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fixedMeanOutletInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    meanValue_->writeData(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
