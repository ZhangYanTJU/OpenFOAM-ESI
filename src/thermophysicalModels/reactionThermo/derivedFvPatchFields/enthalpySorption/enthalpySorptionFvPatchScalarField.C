/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "enthalpySorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "rhoReactionThermo.H"
#include "speciesSorptionFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::enthalpySorptionFvPatchScalarField::enthalpyModelType
>
Foam::enthalpySorptionFvPatchScalarField::enthalpyModelTypeNames
({
    { enthalpyModelType::estimated, "estimated" },
    { enthalpyModelType::calculated, "calculated" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    enthalpyModel_(enthalpyModelType::estimated),
    C_(scalar(0)),
    enthalpyMassLoad_(),
    dhdt_(p.size(), Zero),
    includeHs_(false),
    pName_("p"),
    TName_("T"),
    Hvap_(Zero)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    enthalpyModel_(enthalpyModelTypeNames.get("enthalpyModel", dict)),
    C_(dict.getCheck<scalar>("C", scalarMinMax::ge(0))),
    enthalpyMassLoad_(),
    speciesName_(dict.get<word>("species")),
    dhdt_
    (
        dict.found("dhdt")
      ? scalarField("dhdt", dict, p.size())
      : scalarField(p.size(), 0)
    ),
    includeHs_(dict.getOrDefault<bool>("includeHs", "true")),
    pName_(dict.getOrDefault<word>("p", "p")),
    TName_(dict.getOrDefault<word>("T", "T")),
    Hvap_(dict.getCheck<scalar>("Hvap", scalarMinMax::ge(0)))
{
    switch (enthalpyModel_)
    {
        case enthalpyModelType::calculated:
        {
            enthalpyMassLoad_ = Function1<scalar>::New("enthalpyTable", dict);
            break;
        }
        case enthalpyModelType::estimated:
        {
             break;
        }
    }

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
}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    enthalpyModel_(ptf.enthalpyModel_),
    C_(ptf.C_),
    enthalpyMassLoad_(ptf.enthalpyMassLoad_.clone()),
    speciesName_(ptf.speciesName_),
    dhdt_(ptf.dhdt_),
    includeHs_(ptf.includeHs_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    Hvap_(ptf.Hvap_)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
    enthalpyModel_(ptf.enthalpyModel_),
    C_(ptf.C_),
    enthalpyMassLoad_(ptf.enthalpyMassLoad_.clone()),
    speciesName_(ptf.speciesName_),
    dhdt_(ptf.dhdt_),
    includeHs_(ptf.includeHs_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    Hvap_(ptf.Hvap_)
{}


Foam::enthalpySorptionFvPatchScalarField::enthalpySorptionFvPatchScalarField
(
    const enthalpySorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF),
    C_(ptf.C_),
    enthalpyMassLoad_(ptf.enthalpyMassLoad_.clone()),
    speciesName_(ptf.speciesName_),
    dhdt_(ptf.dhdt_),
    includeHs_(ptf.includeHs_),
    pName_(ptf.pName_),
    TName_(ptf.TName_),
    Hvap_(ptf.Hvap_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enthalpySorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchScalarField::autoMap(m);
    dhdt_.autoMap(m);
}


void Foam::enthalpySorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    const enthalpySorptionFvPatchScalarField& tiptf =
        refCast<const enthalpySorptionFvPatchScalarField>(ptf);

    dhdt_.rmap(tiptf.dhdt_, addr);
}


Foam::tmp<Foam::scalarField> Foam::enthalpySorptionFvPatchScalarField::
patchSource() const
{
    const auto& Yp =
        refCast<const speciesSorptionFvPatchScalarField>
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                speciesName_
            )
        );

    //mass rate [Kg/sec/m3]
    tmp<scalarField> tmassb(Yp.patchSource());
    const scalarField& massb = tmassb();

    // The moles absorbed by the solid
    // dhdt[J/Kg] * Kg/sec/m3 = [J/m3/s]
    scalarField dhdt(dhdt_*massb);

    if (includeHs_)
    {
        const fvPatchField<scalar>& pp =
            patch().lookupPatchField<volScalarField, scalar>(pName_);

        const fvPatchField<scalar>& Tp =
            patch().lookupPatchField<volScalarField, scalar>(TName_);

        const auto& thermo = db().lookupObject<rhoReactionThermo>
        (
            basicThermo::dictName
        );

        const basicSpecieMixture& composition = thermo.composition();

        const label speicesId =
            thermo.composition().species()[speciesName_];

        scalarField hsp(this->patch().size(), Zero);

        forAll(pp, facei)
        {
            hsp[facei] = composition.Hs(speicesId, pp[facei], Tp[facei]);
        }

        dhdt += hsp*massb;
    }

    if (debug)
    {
        Info<< " Patch enthalpy rate min/max [J/m3/sec]: "
            << gMin(dhdt) << " - " << gMax(dhdt) << endl;
    }

    return tmp<scalarField>(new scalarField(dhdt));
}


void Foam::enthalpySorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& Yp =
        refCast<const speciesSorptionFvPatchScalarField>
        (
            patch().lookupPatchField<volScalarField, scalar>
            (
                speciesName_
            )
        );

    switch (enthalpyModel_)
    {
        case enthalpyModelType::estimated:
        {
            dhdt_ = -C_*Hvap_;
            break;
        }
        case enthalpyModelType::calculated:
        {
            // mass [mol/Kg]
            tmp<scalarField> tmassb(Yp.mass());
            const scalarField& massb = tmassb.ref();
            forAll(massb, faceI)
            {
                scalar mFaceI = massb[faceI];

                dhdt_[faceI] = enthalpyMassLoad_->value(mFaceI);
            }
            break;
        }
        default:
            break;
    }

    if (debug)
    {
        Info<< "  Enthalpy change min/max [J/Kg]: "
            << gMin(dhdt_) << " - " << gMax(dhdt_) << endl;
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::enthalpySorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry
    (
        "enthalpyModel", enthalpyModelTypeNames[enthalpyModel_]
    );

    os.writeEntry("C", C_);

    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntryIfDifferent<word>("T", "T", TName_);

    os.writeEntry("Hvap", Hvap_);

    os.writeEntryIfDifferent<bool>("includeHs", true, includeHs_);

    dhdt_.writeEntry("dhdt", os);

    if (enthalpyMassLoad_)
    {
        enthalpyMassLoad_->writeData(os);
    }

    os.writeEntry<word>("species", speciesName_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        enthalpySorptionFvPatchScalarField
    );
}

// ************************************************************************* //
