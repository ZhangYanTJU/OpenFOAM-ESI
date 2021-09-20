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

#include "slowSorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::slowSorptionFvPatchScalarField::slowSorptionModelType
>
Foam::slowSorptionFvPatchScalarField::slowSorptionModelTypeNames
({
    { slowSorptionModelType::HENRY, "Henry" },
    { slowSorptionModelType::LANGMUIR, "Langmuir" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slowSorptionFvPatchScalarField::
slowSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    slowSorptionModel_(slowSorptionModelType::HENRY),
    D_(scalar(1)),
    kAds_(scalar(1)),
    max_(scalar(1)),
    kDes_(scalar(0))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::slowSorptionFvPatchScalarField::
slowSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    slowSorptionModel_(slowSorptionModelTypeNames.get("model", dict)),
    D_(dict.getCheck<scalar>("D", scalarMinMax::ge(SMALL))),
    kAds_(dict.getCheck<scalar>("kAds", scalarMinMax::ge(SMALL))),
    max_(dict.getCheck<scalar>("max", scalarMinMax::ge(SMALL))),
    kDes_(dict.getCheckOrDefault<scalar>("kDes", 0, scalarMinMax::ge(SMALL)))
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}


Foam::slowSorptionFvPatchScalarField::
slowSorptionFvPatchScalarField
(
    const slowSorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    slowSorptionModel_(ptf.slowSorptionModel_),
    D_(ptf.D_),
    kAds_(ptf.kAds_),
    max_(ptf.max_),
    kDes_(ptf.kDes_)
{}


Foam::slowSorptionFvPatchScalarField::
slowSorptionFvPatchScalarField
(
    const slowSorptionFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    slowSorptionModel_(ptf.slowSorptionModel_),
    D_(ptf.D_),
    kAds_(ptf.kAds_),
    max_(ptf.max_),
    kDes_(ptf.kDes_)
{}


Foam::slowSorptionFvPatchScalarField::
slowSorptionFvPatchScalarField
(
    const slowSorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    slowSorptionModel_(ptf.slowSorptionModel_),
    D_(ptf.D_),
    kAds_(ptf.kAds_),
    max_(ptf.max_),
    kDes_(ptf.kDes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slowSorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::slowSorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::slowSorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    switch (slowSorptionModel_)
    {
        case slowSorptionModelType::HENRY:
        {
            // (P:Eq. 2.55) casted into mixed BC form
            valueFraction() =
                (kDes_ - kAds_*max_)
               /(D_*max_*patch().deltaCoeffs() - kDes_ + SMALL);

            break;
        }
        case slowSorptionModelType::LANGMUIR:
        {
            const scalarField co(patchInternalField());

            // (P:Eq. 2.56) casted into mixed BC form
            valueFraction() =
                (kDes_*D_ + kAds_*patch().deltaCoeffs()*(co - max_))
               /(D_*max_ + kDes_*D_ + kAds_*co*patch().deltaCoeffs() + SMALL);

            break;
        }
        default:
            break;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::slowSorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("model", slowSorptionModelTypeNames[slowSorptionModel_]);
    os.writeEntry("D", D_);
    os.writeEntry("kAds", kAds_);
    os.writeEntry("max", max_);
    os.writeEntry("kDes", kDes_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        slowSorptionFvPatchScalarField
    );
}

// ************************************************************************* //
