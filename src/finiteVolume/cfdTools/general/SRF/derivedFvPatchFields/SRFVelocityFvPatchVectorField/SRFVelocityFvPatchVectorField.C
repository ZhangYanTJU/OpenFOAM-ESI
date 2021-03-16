/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "SRFVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "SRFModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    relative_(false),
    inletValue_(p.size(), Zero),
    refValuePtr_(fvPatchVectorField::New("refValue", p, iF))
{}


Foam::SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    inletValue_(ptf.inletValue_, mapper),
    refValuePtr_()
{
    if (ptf.refValuePtr_.valid())
    {
        refValuePtr_ =
            fvPatchVectorField::New(ptf.refValuePtr_(), p, iF, mapper);
    }
}


Foam::SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    relative_(dict.get<Switch>("relative")),
    inletValue_(p.size(), Zero),
    refValuePtr_()
{
    if (dict.found("inletValue"))
    {
        inletValue_ = vectorField("inletValue", dict, p.size());
    }
    else
    {
        if (dict.found("refValue"))
        {
            refValuePtr_ =
                fvPatchVectorField::New(p, iF, dict.subDict("refValue"));
        }
    }
}


Foam::SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& srfvpvf
)
:
    fixedValueFvPatchVectorField(srfvpvf),
    relative_(srfvpvf.relative_),
    inletValue_(srfvpvf.inletValue_),
    refValuePtr_()
{
    if (srfvpvf.refValuePtr_.valid())
    {
        refValuePtr_ = srfvpvf.refValuePtr_().clone();
    }
}


Foam::SRFVelocityFvPatchVectorField::SRFVelocityFvPatchVectorField
(
    const SRFVelocityFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(srfvpvf, iF),
    relative_(srfvpvf.relative_),
    inletValue_(srfvpvf.inletValue_),
    refValuePtr_()
{
    if (srfvpvf.refValuePtr_.valid())
    {
        refValuePtr_ = srfvpvf.refValuePtr_().clone();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SRFVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
    inletValue_.autoMap(m);

    if (refValuePtr_.valid())
    {
        refValuePtr_->autoMap(m);
    }
}


void Foam::SRFVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const auto& srfptf = refCast<const SRFVelocityFvPatchVectorField>(ptf);

    inletValue_.rmap(srfptf.inletValue_, addr);

    if (srfptf.refValuePtr_.valid())
    {
        refValuePtr_->rmap(srfptf.refValuePtr_(), addr);
    }
}


void Foam::SRFVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If not relative to the SRF include the effect of the SRF
    if (!relative_)
    {
        // Get reference to the SRF model
        const auto& srf = db().lookupObject<SRF::SRFModel>("SRFProperties");

        // Determine patch velocity due to SRF
        const vectorField SRFVelocity(srf.velocity(patch().Cf()));

        if (refValuePtr_.valid())
        {
            refValuePtr_->evaluate();
            operator==(-SRFVelocity + refValuePtr_());
        }
        else
        {
            operator==(-SRFVelocity + inletValue_);
        }
    }
    // If already relative to the SRF simply supply the inlet value as a fixed
    // value
    else
    {
        if (refValuePtr_.valid())
        {
            refValuePtr_->evaluate();
            operator==(refValuePtr_());
        }
        else
        {
            operator==(inletValue_);
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::SRFVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("relative", relative_);

    if (refValuePtr_.valid())
    {
        os.beginBlock("refValue");
        refValuePtr_->write(os);
        os.endBlock();
    }
    else
    {
        inletValue_.writeEntry("inletValue", os);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        SRFVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
