/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "transportTwoEqs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(transportTwoEqs, 0);
    addToRunTimeSelectionTable(XiModel, transportTwoEqs, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::transportTwoEqs::transportTwoEqs
(
    const dictionary& XiProperties,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su,
    const volScalarField& rho,
    const volScalarField& b,
    const surfaceScalarField& phi
)
:
    XiModel(XiProperties, thermo, turbulence, Su, rho, b, phi),
    XiShapeCoef_(XiModelCoeffs_.get<scalar>("XiShapeCoef")),
    CpfiDot_(XiModelCoeffs_.get<scalar>("CpfiDot")),
    CpfiCross_(XiModelCoeffs_.get<scalar>("CpfiCross")),
    GEtaExp_(XiModelCoeffs_.get<scalar>("GEtaExp")),
    LOverCw_(XiModelCoeffs_.get<scalar>("LOverCw")),
    XiEqModel_
    (
        XiEqModel::New(XiProperties, "XiEqModel", thermo, turbulence, Su)
    ),
    XiGModel_(XiGModel::New(XiProperties, "XiGModel", thermo, turbulence, Su)),
    XpEqModel_
    (
        XiEqModel::New(XiProperties, "XpEqModel", thermo, turbulence, Su)
    ),
    XpGModel_
    (
        XiGModel::New(XiProperties, "XpGModel", thermo, turbulence, Su)
    ),
    Ep_
    (
        IOobject
        (
            "Ep",
            b.time().timeName(),
            b.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        b.mesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiModels::transportTwoEqs::~transportTwoEqs()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::XiModels::transportTwoEqs::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::transportTwoEqs::correct
(
    const fv::convectionScheme<scalar>& mvConvection
)
{
    const volScalarField XiEqEta(XiEqModel_->XiEq());
    volScalarField GEta(XiGModel_->G());

    GEta *= max( 1.0, exp( GEtaExp_*(1.0 - (Xi_ - 1.0)/(XiEqEta - 0.999)))) ;

    const volScalarField R(GEta*XiEqEta/(XiEqEta - 0.999));

    const volScalarField XiEqStar(R/(R - GEta));

    const volScalarField XiEq
    (
        1.0 + (1.0 + (2*XiShapeCoef_)*(0.5 - b_))*(XiEqStar - 1.0)
    );

    const volScalarField G(R*(XiEq - 1.0)/XiEq);


    const objectRegistry& db = b_.db();
    const volScalarField& betav = db.lookupObject<volScalarField>("betav");
    const volScalarField& p = db.lookupObject<volScalarField>("p");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const surfaceScalarField& phiSt =
        db.lookupObject<surfaceScalarField>("phiSt");
    const volScalarField& Db = db.lookupObject<volScalarField>("Db");
    const surfaceScalarField& nf = db.lookupObject<surfaceScalarField>("nf");

    surfaceScalarField phiXi
    (
        "phiXi",
        phiSt
      + (
          - fvc::interpolate(fvc::laplacian(Db, b_)/mgb)*nf
          + fvc::interpolate(rho_)
          * fvc::interpolate(Su_*(1.0/(Xi_*Xp_) - (Xi_*Xp_)))*nf
        )
    );


    dimensionedScalar zero
    (
        "zero",
        dimensionSet(2, -6, -2, 0, 0, 0, 0),
        scalar(0.0)
    );

    const volScalarField Gpfi
    (
          CpfiDot_*sqrt(max(fvc::grad(rho_)&fvc::grad(p), zero ))/rho_*b_*(1.0-b_)
        + CpfiCross_*sqrt(mag(fvc::grad(rho_)^fvc::grad(p) ))/rho_*b_*(1.0-b_)
    );

    fvScalarMatrix  XiEqn_
    (
        betav*fvm::ddt(rho_, Xi_)
      + mvConvection.fvmDiv(phi_, Xi_)
      + fvm::div(phiXi, Xi_)
      - fvm::Sp(fvc::div(phiXi), Xi_)
     ==
        betav*rho_*(R + Gpfi )
      - fvm::Sp(betav*rho_*(R - G), Xi_)
    );

    XiEqn_.relax();
    XiEqn_.solve();

    // Correct boundedness of Xi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    Xi_.max(1.0);
    Xi_ = min(Xi_, 2.0*XiEq);

    // Calculation of Xp generated by obstacles
    volScalarField XpEqEta("XpEqEta",XpEqModel_->XiEq());

    const volScalarField GpEta("GpEta", XpGModel_->G());

    const volScalarField Rp("Rp", GpEta*XpEqEta/(XpEqEta - 0.999));

     const volScalarField XpEq
    (
        "XpEq",
        1.0 + (1.0 + (2*XiShapeCoef_)*(0.5 - b_))*(XpEqEta - 1.0)
        );

    const volScalarField Gpp("Gpp", Rp*(XpEq - 1.0)/XpEq);


    fvScalarMatrix XpEqn_
    (
        betav*fvm::ddt(rho_, Xp_)
      + mvConvection.fvmDiv(phi_, Xp_)
      + fvm::div(phiXi, Xp_)
      - fvm::Sp(fvc::div(phiXi), Xp_)
     ==
        betav*rho_*Rp
      - fvm::Sp(betav*rho_*(Rp - Gpp), Xp_)
    );

    XpEqn_.relax();
    XpEqn_.solve();

    Xp_.max(1.0);
    Xp_ = min(Xp_, 20.0*XpEq);

    // Calculate Ep
    const volScalarField& Lobs = db.lookupObject<volScalarField>("Lobs");
    const scalarField Cw = pow(Su_.mesh().V(), 2.0/3.0);
    volScalarField LI(Lobs);

    LI.primitiveFieldRef() = max(LI.primitiveField(),LOverCw_*sqrt(Cw));

    fvScalarMatrix EpEqn_
    (
        betav*fvm::ddt(rho_, Ep_)
      + mvConvection.fvmDiv(phi_, Ep_)
      + fvm::div(phiXi, Ep_)
      - fvm::Sp(fvc::div(phiXi), Ep_)
     ==
        betav*rho_*Gpp*Xp_/LI
      - fvm::Sp(betav*rho_*Rp, Ep_)
   );

    EpEqn_.relax();
    EpEqn_.solve();

    Ep_.max(0.0);
    Ep_.min(100000.0);
}


bool Foam::XiModels::transportTwoEqs::read(const dictionary& XiProperties)
{
    XiModel::read(XiProperties);

    XiModelCoeffs_.readEntry("XiShapeCoef", XiShapeCoef_);
    XiModelCoeffs_.readEntry("CpfiDot", CpfiDot_);
    XiModelCoeffs_.readEntry("CpfiCross", CpfiCross_);
    XiModelCoeffs_.readEntry("GEtaExp", GEtaExp_);

    return true;
}


// ************************************************************************* //
