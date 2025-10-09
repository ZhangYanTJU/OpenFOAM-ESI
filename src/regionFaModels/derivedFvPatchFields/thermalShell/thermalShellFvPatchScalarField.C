/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
#include "dictionaryContent.H"
#include "regionProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void thermalShellFvPatchScalarField::create_baffle()
{
    if (!baffle_)
    {
        baffle_.reset
        (
            baffleType::New(this->patch().boundaryMesh().mesh(), dict_)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
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
    parent_bctype(ptf, p, iF, mapper),
    dict_(ptf.dict_)
{}


thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
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
    // Create baffle now.
    // Lazy evaluation has issues with loading the finite-area fields
    if (regionModels::allowFaModels())
    {
        create_baffle();
    }
}


thermalShellFvPatchScalarField::thermalShellFvPatchScalarField
(
    const thermalShellFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    dict_(ptf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalShellFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField& pfld = *this;

    // Create baffle if needed
    if (!baffle_)
    {
        create_baffle();
    }

    auto& baffle = baffle_();

    baffle.evolve();

    baffle.vsm().mapToVolumePatch(baffle.T(), pfld, patch().index());

    parent_bctype::updateCoeffs();
}


void thermalShellFvPatchScalarField::write(Ostream& os) const
{
    parent_bctype::write(os);
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
