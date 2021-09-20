/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "fastSorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fastSorptionFvPatchScalarField::fastSorptionModelType
>
Foam::fastSorptionFvPatchScalarField::fastSorptionModelTypeNames
({
    { fastSorptionModelType::HENRY, "Henry" },
    { fastSorptionModelType::LANGMUIR, "Langmuir" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fastSorptionFvPatchScalarField::
fastSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fastSorptionModel_(fastSorptionModelType::HENRY),
    k_(scalar(1)),
    max_(scalar(1)),
    aL_(scalar(1))
{}


Foam::fastSorptionFvPatchScalarField::
fastSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    fastSorptionModel_(fastSorptionModelTypeNames.get("model", dict)),
    k_(scalar(1)),
    max_(scalar(1)),
    aL_(scalar(1))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(Zero);
    }

    const word& modelType = fastSorptionModelTypeNames[fastSorptionModel_];

    const dictionary& coeffs = dict.optionalSubDict(modelType + "Coeffs");

    switch (fastSorptionModel_)
    {
        case fastSorptionModelType::HENRY:
        {
            k_ = coeffs.getCheck<scalar>("k", scalarMinMax::ge(SMALL));
            break;
        }
        case fastSorptionModelType::LANGMUIR:
        {
            max_ = coeffs.getCheck<scalar>("max", scalarMinMax::ge(SMALL));
            aL_ = coeffs.getCheck<scalar>("aL", scalarMinMax::ge(SMALL));
            break;
        }
        default:
            break;
    }
}


Foam::fastSorptionFvPatchScalarField::
fastSorptionFvPatchScalarField
(
    const fastSorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fastSorptionModel_(ptf.fastSorptionModel_),
    k_(ptf.k_),
    max_(ptf.max_),
    aL_(ptf.aL_)
{}


Foam::fastSorptionFvPatchScalarField::
fastSorptionFvPatchScalarField
(
    const fastSorptionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    fastSorptionModel_(ptf.fastSorptionModel_),
    k_(ptf.k_),
    max_(ptf.max_),
    aL_(ptf.aL_)
{}


Foam::fastSorptionFvPatchScalarField::
fastSorptionFvPatchScalarField
(
    const fastSorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    fastSorptionModel_(ptf.fastSorptionModel_),
    k_(ptf.k_),
    max_(ptf.max_),
    aL_(ptf.aL_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fastSorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::fastSorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::fastSorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField co(patchInternalField());

    switch (fastSorptionModel_)
    {
        case fastSorptionModelType::HENRY:
        {
            // (P:Eq. 2.39)
            operator==(k_*co);
            break;
        }
        case fastSorptionModelType::LANGMUIR:
        {
            // (P:Eq. 2.41)
            operator==(max_*co/(aL_ + co));
            break;
        }
        default:
            break;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fastSorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    const word& modelType = fastSorptionModelTypeNames[fastSorptionModel_];

    os.writeEntry("model", modelType);

    os.beginBlock(word(modelType + "Coeffs"));

    switch (fastSorptionModel_)
    {
        case fastSorptionModelType::HENRY:
        {
            os.writeEntry("k", k_);
            break;
        }
        case fastSorptionModelType::LANGMUIR:
        {
            os.writeEntry("max", max_);
            os.writeEntry("aL", aL_);
            break;
        }
        default:
            break;
    }

    os.endBlock();

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fastSorptionFvPatchScalarField
    );
}

// ************************************************************************* //
