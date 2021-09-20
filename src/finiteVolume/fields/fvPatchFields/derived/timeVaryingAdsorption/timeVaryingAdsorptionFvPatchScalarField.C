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

#include "timeVaryingAdsorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::timeVaryingAdsorptionFvPatchScalarField::ddtSchemeType
>
Foam::timeVaryingAdsorptionFvPatchScalarField::ddtSchemeTypeNames_
({
    {
        ddtSchemeType::tsEuler,
        fv::EulerDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsCrankNicolson,
        fv::CrankNicolsonDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsBackward,
        fv::backwardDdtScheme<scalar>::typeName_()
    },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingAdsorptionFvPatchScalarField::
timeVaryingAdsorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    k_(scalar(1)),
    max_(scalar(1))
{}


Foam::timeVaryingAdsorptionFvPatchScalarField::
timeVaryingAdsorptionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    k_(dict.getCheck<scalar>("k", scalarMinMax::ge(SMALL))),
    max_(dict.getCheck<scalar>("max", scalarMinMax::ge(SMALL)))
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


Foam::timeVaryingAdsorptionFvPatchScalarField::
timeVaryingAdsorptionFvPatchScalarField
(
    const timeVaryingAdsorptionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    k_(ptf.k_),
    max_(ptf.max_)
{}


Foam::timeVaryingAdsorptionFvPatchScalarField::
timeVaryingAdsorptionFvPatchScalarField
(
    const timeVaryingAdsorptionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    k_(ptf.k_),
    max_(ptf.max_)
{}


Foam::timeVaryingAdsorptionFvPatchScalarField::
timeVaryingAdsorptionFvPatchScalarField
(
    const timeVaryingAdsorptionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    k_(ptf.k_),
    max_(ptf.max_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingAdsorptionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::timeVaryingAdsorptionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::timeVaryingAdsorptionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const scalar dt = db().time().deltaTValue();

    const auto& fld =
        db().lookupObject<volScalarField>(this->internalField().name());
    const volScalarField& fld0 = fld.oldTime();

    // Lookup d/dt scheme from database
    const word ddtSchemeName(fld.mesh().ddtScheme(fld.name()));
    const ddtSchemeType& ddtScheme = ddtSchemeTypeNames_[ddtSchemeName];


    tmp<scalarField> dfldp = k_*(max_ - patchInternalField())*dt;

    switch (ddtScheme)
    {
        case tsEuler:
        case tsCrankNicolson:
        {
            operator==(fld0.boundaryField()[patchi] + dfldp);

            break;
        }
        case tsBackward:
        {
            const scalar dt0 = db().time().deltaT0Value();

            const scalar c = scalar(1) + dt/(dt + dt0);
            const scalar c00 = dt*dt/(dt0*(dt + dt0));
            const scalar c0 = c + c00;

            operator==
                (
                    (
                        c0*fld0.boundaryField()[patchi]
                    - c00*fld0.oldTime().boundaryField()[patchi]
                    + dfldp
                    )/c
                );

            break;
        }
        default:
        {
            FatalErrorInFunction
                << ddtSchemeName << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingAdsorptionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("k", k_);
    os.writeEntry("max", max_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingAdsorptionFvPatchScalarField
    );
}

// ************************************************************************* //
