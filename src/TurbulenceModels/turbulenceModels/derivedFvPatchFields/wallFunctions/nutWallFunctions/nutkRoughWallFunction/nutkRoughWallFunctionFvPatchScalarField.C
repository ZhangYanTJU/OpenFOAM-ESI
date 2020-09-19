/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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

#include "nutkRoughWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::nutkRoughWallFunctionFvPatchScalarField::fnRough
(
    const scalar KsPlus,
    const scalar Cs
) const
{
    // Return fn based on non-dimensional roughness height
    if (KsPlus < 90.0)
    {
        return pow
        (
            (KsPlus - 2.25)/87.75 + Cs*KsPlus,
            sin(0.4258*(log(KsPlus) - 0.811))
        );
    }

    return (1.0 + Cs*KsPlus);
}


Foam::tmp<Foam::scalarField> Foam::nutkRoughWallFunctionFvPatchScalarField::
calcNut() const
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

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const auto& nuw = tnuw();

    uStar_ = calcUStar();
    yPlus_ = calcYPlus();

    auto tnutw = tmp<scalarField>::New(*this);
    auto& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        const scalar KsPlus = uStar_[facei]*Ks_[facei]/nuw[facei];

        scalar Edash = E_;
        if (2.25 < KsPlus)
        {
            Edash /= fnRough(KsPlus, Cs_[facei]);
        }

        const scalar limitingNutw = max(nutw[facei], nuw[facei]);

        // To avoid oscillations limit the change in the wall viscosity
        // which is particularly important if it temporarily becomes zero
        nutw[facei] =
            max
            (
                min
                (
                    nuw[facei]
                   *(
                       yPlus_[facei]*kappa_
                      /log(max(Edash*yPlus_[facei], 1+1e-4)) - 1
                    ),
                    2*limitingNutw
                ), 0.5*limitingNutw
            );

        if (debug)
        {
            Info<< "yPlus = " << yPlus_[facei]
                << ", KsPlus = " << KsPlus
                << ", Edash = " << Edash
                << ", nutw = " << nutw[facei]
                << endl;
        }
    }

    return tnutw;
}


Foam::scalarField Foam::nutkRoughWallFunctionFvPatchScalarField::
calcUStar() const
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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();

    return pow025(Cmu_)*sqrt(kwc);
}


Foam::scalarField Foam::nutkRoughWallFunctionFvPatchScalarField::
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

    return uStar_*turbModel.y()[patchi]/turbModel.nu(patchi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    Ks_(p.size(), Zero),
    Cs_(p.size(), Zero),
    uStar_(p.size(), Zero)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ks_(ptf.Ks_, mapper),
    Cs_(ptf.Cs_, mapper),
    uStar_(ptf.uStar_, mapper)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    Cs_("Cs", dict, p.size()),
    uStar_(p.size(), Zero)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_),
    uStar_(rwfpsf.uStar_)
{}


Foam::nutkRoughWallFunctionFvPatchScalarField::
nutkRoughWallFunctionFvPatchScalarField
(
    const nutkRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(rwfpsf, iF),
    Ks_(rwfpsf.Ks_),
    Cs_(rwfpsf.Cs_),
    uStar_(rwfpsf.uStar_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nutkRoughWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    nutkWallFunctionFvPatchScalarField::autoMap(m);
    Ks_.autoMap(m);
    Cs_.autoMap(m);
}


void Foam::nutkRoughWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutkWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutkRoughWallFunctionFvPatchScalarField& nrwfpsf =
        refCast<const nutkRoughWallFunctionFvPatchScalarField>(ptf);

    Ks_.rmap(nrwfpsf.Ks_, addr);
    Cs_.rmap(nrwfpsf.Cs_, addr);
}


void Foam::nutkRoughWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    Cs_.writeEntry("Cs", os);
    Ks_.writeEntry("Ks", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutkRoughWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
