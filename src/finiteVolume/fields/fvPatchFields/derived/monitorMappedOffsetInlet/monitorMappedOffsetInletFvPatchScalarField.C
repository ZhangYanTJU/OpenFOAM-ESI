/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "monitorMappedOffsetInletFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monitorMappedOffsetInletFvPatchScalarField::
monitorMappedOffsetInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    monitorPatchName_(),
    phiName_("phi"),
    offset_(0),
    max_(GREAT),
    min_(SMALL)
{}


Foam::monitorMappedOffsetInletFvPatchScalarField::
monitorMappedOffsetInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    monitorPatchName_(dict.get<word>("monitorPatch")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    offset_(dict.get<scalar>("offset")),
    max_(dict.getOrDefault<scalar>("max", GREAT)),
    min_(dict.getOrDefault<scalar>("min", SMALL))
{
    if (max_ < min_)
    {
        FatalIOErrorInFunction(dict)
            << "Upper-limit scalar max = " << max_ << nl
            << " cannot be smaller than lower-limit scalar min = " << min_
            << exit(FatalIOError);
    }

    if (max_ < 0 || min_ < 0)
    {
        WarningInFunction
            << "max = " << max_ << ", or min = " << min_ << " is negative"
            << endl;
    }
}


Foam::monitorMappedOffsetInletFvPatchScalarField::
monitorMappedOffsetInletFvPatchScalarField
(
    const monitorMappedOffsetInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    offset_(ptf.offset_),
    max_(ptf.max_),
    min_(ptf.min_)
{}


Foam::monitorMappedOffsetInletFvPatchScalarField::
monitorMappedOffsetInletFvPatchScalarField
(
    const monitorMappedOffsetInletFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    offset_(ptf.offset_),
    max_(ptf.max_),
    min_(ptf.min_)
{}


Foam::monitorMappedOffsetInletFvPatchScalarField::
monitorMappedOffsetInletFvPatchScalarField
(
    const monitorMappedOffsetInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    offset_(ptf.offset_),
    max_(ptf.max_),
    min_(ptf.min_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::monitorMappedOffsetInletFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const auto& f = dynamic_cast<const volScalarField&>(this->internalField());

    const fvPatch& p = this->patch();
    const label monitorPatchID =
        p.patch().boundaryMesh().findPatchID(monitorPatchName_);

    if (monitorPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find monitor patch " << monitorPatchName_
            << exit(FatalError);
    }

    const fvPatch& monitorPatch = p.boundaryMesh()[monitorPatchID];

    const fvPatchScalarField& monitorPatchField =
        f.boundaryField()[monitorPatchID];

    const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& monitorPatchPhi = phi.boundaryField()[monitorPatchID];
    const scalar sumMonitorPatchPhi = gSum(monitorPatchPhi);

    if (sumMonitorPatchPhi > SMALL)
    {
        const scalar averageMonitorField =
            gSum(monitorPatchPhi*monitorPatchField)/sumMonitorPatchPhi;

        operator==(min(max(averageMonitorField + offset_, min_), max_));
    }
    else
    {
        const scalar averageMonitorField =
            gSum(monitorPatch.magSf()*monitorPatchField)
           /gSum(monitorPatch.magSf());

        operator==(averageMonitorField);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::monitorMappedOffsetInletFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("monitorPatch", monitorPatchName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("offset", offset_);
    os.writeEntry("max", max_);
    os.writeEntry("min", min_);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        monitorMappedOffsetInletFvPatchScalarField
    );
}


// ************************************************************************* //
