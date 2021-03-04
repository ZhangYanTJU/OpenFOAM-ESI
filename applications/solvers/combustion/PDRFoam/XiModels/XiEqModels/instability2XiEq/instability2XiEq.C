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

#include "instability2XiEq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(instability2XiEq, 0);
    addToRunTimeSelectionTable(XiEqModel, instability2XiEq, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::instability2XiEq::instability2XiEq
(
    const dictionary& XiEqProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, modelType, thermo, turbulence, Su),
    saModel_
    (
        IOdictionary
        (
            IOobject
            (
                "combustionProperties",
                Su.mesh().time().constant(),
                Su.mesh(),
                IOobject::MUST_READ
            )
        ),
        thermo
    ),
    CIn_(saModel_.CIn()),
    defaultCIn_(XiEqModelCoeffs_.get<scalar>("defaultCIn")),
    XiEqInFade_(XiEqModelCoeffs_.get<scalar>("XiEqInFade")),
    XiEqModel_
    (
        XiEqModel::New(XiEqModelCoeffs_, modelType, thermo, turbulence, Su)
    )
{
    if (CIn_ <= 0.0)
    {
        CIn_ = defaultCIn_;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::instability2XiEq::XiEq() const
{
    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            Su_.mesh().time().constant(),
            Su_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    ignition ign(combustionProperties, Su_.mesh().time(), Su_.mesh());

    //const scalar ignTim = ign.sites()[0].tmIgn();
    const scalar curTime = Su_.mesh().time().value();
    const scalar deltaT = Su_.mesh().time().deltaTValue();
    const scalar ignTim = curTime - deltaT - ign.sites()[0].time();

    volScalarField turbXiEq(XiEqModel_->XiEq());

    volScalarField XiEqIn1("XiEqIn1", 0.0*turbXiEq);

    dimensionedScalar CIn("CIn", dimensionSet(0, -2, 1, 0, 0, 0, 0), CIn_);
    dimensionedScalar ignTm("ignTm", dimTime, ignTim);
    XiEqIn1 =  exp(CIn*Su_*Su_*ignTm) - 1.0;

    return
    (
        1.0 + sqrt(XiEqInFade_*sqr(XiEqIn1) + sqr(turbXiEq - 1.0))
    );
}


bool Foam::XiEqModels::instability2XiEq::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    XiEqModelCoeffs_.readEntry("defaultCIn", defaultCIn_);
    XiEqModelCoeffs_.readEntry("XiEqInFade", XiEqInFade_);

    return XiEqModel_->read(XiEqModelCoeffs_);
}


// ************************************************************************* //
