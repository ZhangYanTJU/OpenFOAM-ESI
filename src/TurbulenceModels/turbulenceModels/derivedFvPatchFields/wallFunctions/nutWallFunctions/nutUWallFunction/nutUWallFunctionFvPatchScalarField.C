/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "nutUWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    yPlus_ = calcYPlus();

    const scalarField nutLog
    (
        turbModel.nu(patchi)
      * (yPlus_*kappa_/log(max(E_*yPlus_, 1 + 1e-4)) - 1)
    );

    const scalar nutVis = 0;

    return blend(nutLog, nutVis);
}


Foam::scalarField Foam::nutUWallFunctionFvPatchScalarField::
calcYPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const auto& nuw = tnuw();

    // Viscous sublayer estimation
    scalarField yplus(y*sqrt(turbModel.nuEff(patchi)*mag(Uw.snGrad()))/nuw);

    // Inertial sublayer estimation
    scalarField yplusLog(patch().size());
    forAll(yplusLog, facei)
    {
        const scalar kappaRe = kappa_*magUp[facei]*y[facei]/nuw[facei];

        scalar yp = yPlusLam_;
        const scalar ryPlusLam = 1/yp;

        int iter = 0;
        scalar yPlusLast = 0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1 + log(E_*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );

        yplusLog[facei] = max(0, yp);
    }

    // Replace viscous estimation with inertial est. if yplusLog > yPlusLam
    std::copy_if
    (
        yplusLog.cbegin(),
        yplusLog.cend(),
        yplus.begin(),
        [&](const scalar& ypLog){ return ypLog > yPlusLam_; }
    );

    return yplus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nutUWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutUWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
