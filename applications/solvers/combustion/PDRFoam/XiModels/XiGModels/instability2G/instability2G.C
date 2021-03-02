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

#include "IFstream.H"
#include "instability2G.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(instability2G, 0);
    addToRunTimeSelectionTable(XiGModel, instability2G, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::instability2G::instability2G
(
    const dictionary& XiGProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModel(XiGProperties, modelType, thermo, turbulence, Su),
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
    defaultCIn_(XiGModelCoeffs_.get<scalar>("defaultCIn")),
    GInFade_(XiGModelCoeffs_.get<scalar>("GInFade")),
    GInMult_(XiGModelCoeffs_.get<scalar>("GInMult")),
    lambdaIn_("lambdaIn", XiGModelCoeffs_),
    XiGModel_
    (
        XiGModel::New(XiGModelCoeffs_,modelType,thermo, turbulence, Su)
    )
{
    if (CIn_ <=  0.0)
    {
        CIn_ = defaultCIn_;
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiGModels::instability2G::~instability2G()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::instability2G::G() const
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

    scalar curTime = Su_.mesh().time().value();
    scalar deltaT = Su_.mesh().time().deltaTValue();
    const scalar ignTim = curTime - deltaT - ign.sites()[0].time();

    volScalarField turbXiG(XiGModel_->G());

    volScalarField GIn("GIn", 0.0*turbXiG);


    forAll (GIn, i)
    {
        GIn[i] = CIn_*Su_[i]*Su_[i]*exp(CIn_*Su_[i]*Su_[i]*ignTim)*GInMult_;
    }

    dimensionedScalar CIn("CIn", dimensionSet(0, -2, 1, 0, 0, 0, 0), CIn_);
    dimensionedScalar ignTm("ignTm", dimTime, ignTim);

    GIn = CIn*Su_*Su_*exp(CIn*Su_*Su_*ignTm)*GInMult_;

    GIn *=
        GIn
        /
        (
            GIn
          + GInFade_*turbXiG
          + dimensionedScalar("GSmall", inv(dimTime), SMALL)
        );

    return (GIn + turbXiG);
}


Foam::tmp<Foam::volScalarField> Foam::XiGModels::instability2G::Db() const
{
    const objectRegistry& db = Su_.db();
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");
    const volScalarField& rho = db.lookupObject<volScalarField>("rho");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const volScalarField& Db1 = db.lookupObject<volScalarField>("Db");

    //return  turbulence_.muEff()
    return Db1
        + rho*Su_*(Xi - 1.0)*mgb*(0.5*lambdaIn_)/(mgb + 1.0/lambdaIn_);
}


bool Foam::XiGModels::instability2G::read(const dictionary& XiGProperties)
{
    XiGModel::read(XiGProperties);

    XiGModelCoeffs_.readEntry("defaultCIn", defaultCIn_);
    XiGModelCoeffs_.readEntry("GInFade", GInFade_);
    XiGModelCoeffs_.readEntry("GInMult", GInMult_);
    XiGModelCoeffs_.readEntry("lambdaIn", lambdaIn_);

    return true;
}


// ************************************************************************* //
