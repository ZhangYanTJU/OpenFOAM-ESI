/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "vibrationShellFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    baffle_(),
    dict_(dictionary::null)
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const vibrationShellFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>
    (
        ptf,
        p,
        iF,
        mapper
    ),
    baffle_(),
    dict_(ptf.dict_)
{}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    baffle_(),
    dict_(dict)
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
    
    if (baffle_.empty())
    {
        baffle_.reset(regionModels::vibrationShellModel::New(p, dict).ptr());
    }
}


vibrationShellFvPatchScalarField::vibrationShellFvPatchScalarField
(
    const vibrationShellFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    baffle_(),
    dict_(ptf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void vibrationShellFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    baffle_->evolve();

    /*
    volScalarField::Boundary& vfb =
        db().lookupObjectRef<volScalarField>
        (
            this->internalField().name()
        ).boundaryFieldRef();

    baffle_->vsm().mapToVolume(baffle_->T(), vfb);
    */
    
    const areaScalarField aRho(-baffle_->solid().rho()*baffle_->a());

    baffle_->vsm().mapToField(aRho, this->refGrad());

    this->refValue() = Zero;
    this->valueFraction() = Zero;

    mixedFvPatchField<scalar>::updateCoeffs();
}


void vibrationShellFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);

    // Remove value and type already written by fixedValueFvPatchField
    //dict_.remove("value");
    //dict_.remove("type");
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    vibrationShellFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
