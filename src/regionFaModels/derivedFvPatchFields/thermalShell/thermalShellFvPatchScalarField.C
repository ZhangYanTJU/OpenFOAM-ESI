/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "thermalShellFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    bafflePtr_(nullptr),
    dict_()
{}


thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const thermalShellFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField
    (
        ptf,
        p,
        iF,
        mapper
    ),
    bafflePtr_(nullptr),
    dict_(ptf.dict_)
{}


thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    bafflePtr_(nullptr),
    dict_
    (
        // Copy dictionary, but without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value"
            })
        )
    )
{
    if (!bafflePtr_)
    {
        bafflePtr_.reset(baffleType::New(p.boundaryMesh().mesh(), dict_));
    }
}


thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const thermalShellFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    bafflePtr_(nullptr),
    dict_(ptf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalShellFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    bafflePtr_->evolve();

    scalarField& pfld = *this;

    bafflePtr_->vsm().mapToVolumePatch(bafflePtr_->T(), pfld, patch().index());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void thermalShellFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalShellFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
