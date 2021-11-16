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

#include "speciesSorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
const Foam::Enum
<
    Foam::speciesSorptionFvPatchScalarField::equilibriumModelType
>
Foam::speciesSorptionFvPatchScalarField::equilibriumModelTypeNames
({
    { equilibriumModelType::LANGMUIR, "Langmuir" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::speciesSorptionFvPatchScalarField::calcMoleFractions()
{
    tmp<scalarField> tMole(new scalarField(patch().size(), 0));
    scalarField& Mole = tMole.ref();

    if (db().foundObject<rhoReactionThermo>(basicThermo::dictName))
    {
        const auto& thermo =
            db().lookupObject<rhoReactionThermo>
            (
                basicThermo::dictName
            );

        const PtrList<volScalarField>& Y = thermo.composition().Y();

        const volScalarField W(thermo.W());

        const labelUList& faceCells = patch().faceCells();

        const label speicesId =
            thermo.composition().species()[this->internalField().name()];

        const dimensionedScalar Wi
        (
            dimMass/dimMoles,
            thermo.composition().W(speicesId)
        );

        const volScalarField X(W*Y[speicesId]/Wi);

        forAll(faceCells, i)
        {
            label cellId = faceCells[i];
            Mole[i] = X[cellId];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Thermo type is not 'rhoReactionThermo'. " << nl
            << "This BC is designed to operate with a rho based thermo."
            << exit(FatalError);
    }

    return tMole;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(p, iF),
    equilibriumModel_(equilibriumModelType::LANGMUIR),
    kabs_(scalar(1)),
    kl_(0),
    max_(scalar(1)),
    dfldp_(p.size(), 0),
    mass_(p.size(), 0)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    equilibriumModel_(equilibriumModelTypeNames.get("model", dict)),
    kabs_(dict.getCheck<scalar>("kabs", scalarMinMax::ge(0))),
    kl_(dict.getCheck<scalar>("kl", scalarMinMax::ge(0))),
    max_(dict.getCheck<scalar>("max", scalarMinMax::ge(0))),
    dfldp_
    (
        dict.found("dfldp")
      ? scalarField("dfldp", dict, p.size())
      : scalarField(p.size(), 0)
    ),
    mass_
    (
        dict.found("mass")
      ? scalarField("mass", dict, p.size())
      : scalarField(p.size(), 0)
    )
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
}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    equilibriumModel_(ptf.equilibriumModel_),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    dfldp_(ptf.dfldp_, mapper),
    mass_(ptf.mass_, mapper)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf
)
:
    zeroGradientFvPatchScalarField(ptf),
    equilibriumModel_(ptf.equilibriumModel_),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    dfldp_(ptf.dfldp_),
    mass_(ptf.mass_)
{}


Foam::speciesSorptionFvPatchScalarField::speciesSorptionFvPatchScalarField
(
    const speciesSorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(ptf, iF),
    equilibriumModel_(ptf.equilibriumModel_),
    kabs_(ptf.kabs_),
    kl_(ptf.kl_),
    max_(ptf.max_),
    dfldp_(ptf.dfldp_),
    mass_(ptf.mass_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesSorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    zeroGradientFvPatchScalarField::autoMap(m);
    dfldp_.autoMap(m);
    mass_.autoMap(m);
}


void Foam::speciesSorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    zeroGradientFvPatchScalarField::rmap(ptf, addr);

    const speciesSorptionFvPatchScalarField& tiptf =
        refCast<const speciesSorptionFvPatchScalarField>(ptf);

    dfldp_.rmap(tiptf.dfldp_, addr);
    mass_.rmap(tiptf.mass_, addr);
}


Foam::tmp<Foam::scalarField> Foam::speciesSorptionFvPatchScalarField::
patchSource() const
{
    const auto& thermo = db().lookupObject<rhoReactionThermo>
    (
        basicThermo::dictName
    );

    const label speicesId =
        thermo.composition().species()[this->internalField().name()];

    const scalar Wi(thermo.composition().W(speicesId));
    // [mol/Kg/sec]*[g/mol]
    return tmp<scalarField>(new scalarField(-dfldp_*Wi*1e-3));
}


void Foam::speciesSorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar dt = db().time().deltaTValue();

    // mole fraction
    tmp<scalarField> co = calcMoleFractions();

    // equilibrium in mol/Kg
    scalarField cEq(patch().size(), 0);

    switch (equilibriumModel_)
    {
        case equilibriumModelType::LANGMUIR:
        {
            cEq = max_*(kl_*co()/(1 + kl_*co()));
            break;
        }
        default:
            break;
    }

    // source [mol/Kg/sec]
    dfldp_ = kabs_*(max(cEq - mass_, 0.0));

    mass_ += dfldp_*dt;

    if (debug)
    {
        Info<< "  Absorption rate: [mol/Kg/sec] "
            << gMin(dfldp_) << " - " << gMax(dfldp_) << endl;

        Info<< " Absorbed mass max/min: "
            << gMin(mass_) << " - " << gMax(mass_) << endl;
    }

    zeroGradientFvPatchScalarField::updateCoeffs();
}


void Foam::speciesSorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("model", equilibriumModelTypeNames[equilibriumModel_]);
    os.writeEntry("kabs", kabs_);
    os.writeEntry("kl", kl_);
    os.writeEntry("max", max_);
    dfldp_.writeEntry("dfldp", os);
    mass_.writeEntry("mass", os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        speciesSorptionFvPatchScalarField
    );
}

// ************************************************************************* //
