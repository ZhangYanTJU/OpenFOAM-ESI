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

#include "monitorMappedFilterInletFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monitorMappedFilterInletFvPatchScalarField::
monitorMappedFilterInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    monitorPatchName_(),
    phiName_("phi"),
    filterRatio_(0),
    shareRatio_(1)
{}


Foam::monitorMappedFilterInletFvPatchScalarField::
monitorMappedFilterInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    monitorPatchName_(dict.get<word>("monitorPatch")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    filterRatio_(dict.getOrDefault<scalar>("filterRatio", 0)),
    shareRatio_(dict.getOrDefault<scalar>("shareRatio", 1))
{
    if
    (
        max(shareRatio_, filterRatio_) >= 1
     || min(shareRatio_, filterRatio_) < 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "Filtering ratio =" << filterRatio_ << nl
            << "or share ratio =" << shareRatio_ << nl
            << "cannot be larger than 1 or negative"
            << exit(FatalIOError);
    }
}


Foam::monitorMappedFilterInletFvPatchScalarField::
monitorMappedFilterInletFvPatchScalarField
(
    const monitorMappedFilterInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    filterRatio_(ptf.filterRatio_),
    shareRatio_(ptf.shareRatio_)
{}


Foam::monitorMappedFilterInletFvPatchScalarField::
monitorMappedFilterInletFvPatchScalarField
(
    const monitorMappedFilterInletFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    filterRatio_(ptf.filterRatio_),
    shareRatio_(ptf.shareRatio_)
{}


Foam::monitorMappedFilterInletFvPatchScalarField::
monitorMappedFilterInletFvPatchScalarField
(
    const monitorMappedFilterInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    monitorPatchName_(ptf.monitorPatchName_),
    phiName_(ptf.phiName_),
    filterRatio_(ptf.filterRatio_),
    shareRatio_(ptf.shareRatio_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::monitorMappedFilterInletFvPatchScalarField::updateCoeffs()
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
            << abort(FatalError);
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

        operator==(averageMonitorField*(1 - filterRatio_)*shareRatio_);
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


void Foam::monitorMappedFilterInletFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("monitorPatch", monitorPatchName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry("filterRatio", filterRatio_);
    os.writeEntry("shareRatio", shareRatio_);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        monitorMappedFilterInletFvPatchScalarField
    );
}


// ************************************************************************* //
