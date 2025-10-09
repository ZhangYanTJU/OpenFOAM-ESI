/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "velocityFilmShellFvPatchVectorField.H"
#include "dictionaryContent.H"
#include "regionProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void velocityFilmShellFvPatchVectorField::create_baffle()
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

velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF),
    dict_(),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
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
                "value", "refValue", "refGradient", "valueFraction"
            })
        )
    ),
    curTimeIndex_(-1),
    zeroWallVelocity_(dict.getOrDefault<bool>("zeroWallVelocity", true))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1;
    }

    // Create baffle now.
    // Lazy evaluation has issues with loading the finite-area fields
    if (regionModels::allowFaModels())
    {
        create_baffle();
    }
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityFilmShellFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Create baffle if needed
    if (!baffle_)
    {
        create_baffle();
    }

    // Execute the change only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        vectorField& pfld = *this;

        auto& baffle = baffle_();

        baffle.evolve();

        baffle.vsm().mapToVolumePatch(baffle.Us(), pfld, patch().index());

        refGrad() = Zero;
        valueFraction() = 1;

        if (zeroWallVelocity_)
        {
            refValue() = Zero;
        }
        else
        {
            refValue() = pfld;
        }

        curTimeIndex_ = this->db().time().timeIndex();
    }

    parent_bctype::updateCoeffs();
}


void velocityFilmShellFvPatchVectorField::write(Ostream& os) const
{
    parent_bctype::write(os);
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    velocityFilmShellFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
