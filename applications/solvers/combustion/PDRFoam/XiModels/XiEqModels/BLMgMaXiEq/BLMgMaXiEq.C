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

#include "BLMgMaXiEq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(BLMgMaXiEq, 0);
    addToRunTimeSelectionTable(XiEqModel, BLMgMaXiEq, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::BLMgMaXiEq::BLMgMaXiEq
(
    const dictionary& XiEqProperties,
    const word& modelType,
    const psiuReactionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, modelType, thermo, turbulence, Su),
    kaCoef_(XiEqModelCoeffs_.get<scalar>("kaCoef")),
    lowK0_(XiEqModelCoeffs_.get<scalar>("lowK0")),
    lowKg_(XiEqModelCoeffs_.get<scalar>("lowKg")),
    XiEqCoef_(XiEqModelCoeffs_.get<scalar>("XiEqCoef")),
    alphaCoefP_(XiEqModelCoeffs_.get<scalar>("alphaCoefP")),
    betaCoefP_(XiEqModelCoeffs_.get<scalar>("betaCoefP")),
    alphaCoefN_(XiEqModelCoeffs_.get<scalar>("alphaCoefN")),
    betaCoefN_(XiEqModelCoeffs_.get<scalar>("betaCoefN")),
    maLim_(XiEqModelCoeffs_.get<scalar>("maLim")),
    maLim1_(XiEqModelCoeffs_.get<scalar>("maLim1")),
    quenchCoef_(XiEqModelCoeffs_.get<scalar>("quenchCoef")),
    quenchExp_(XiEqModelCoeffs_.get<scalar>("quenchExp")),
    quenchM_(XiEqModelCoeffs_.get<scalar>("quenchM")),
    quenchRate1_(XiEqModelCoeffs_.get<scalar>("quenchRate1")),
    quenchRate2_(XiEqModelCoeffs_.get<scalar>("quenchRate2")),
    lCoef_(XiEqModelCoeffs_.get<scalar>("lCoef")),
    SuMin_(0.01*Su.average()),
    uPrimeCoef_(XiEqModelCoeffs_.get<scalar>("uPrimeCoef")),
    nrExp_(XiEqModelCoeffs_.get<scalar>("nrExp")),
    subGridSchelkin_(XiEqModelCoeffs_.get<bool>("subGridSchelkin")),
    MaModel
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
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModels::BLMgMaXiEq::~BLMgMaXiEq()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::BLMgMaXiEq::XiEq() const
{
    const volScalarField& k = turbulence_.k();
    const volScalarField& epsilon = turbulence_.epsilon();
    volScalarField up("up", sqrt((2.0/3.0)*k));
    if (subGridSchelkin_)
    {
        up.primitiveFieldRef() +=
            calculateSchelkinEffect(uPrimeCoef_, nrExp_);
    }

    volScalarField l(lCoef_*sqrt(3.0/2.0)*up*k/epsilon);
    volScalarField Rl(up*l*thermo_.rhou()/thermo_.muu());

    volScalarField upBySu("upBySu", up/(Su_ + SuMin_));
    volScalarField K("K", kaCoef_*upBySu*upBySu/sqrt(Rl));
    volScalarField Ma("Ma", MaModel.Ma());

    volScalarField regime("regime", MaModel.Ma()*scalar(0.0));

    tmp<volScalarField> tXiEq
    (
        new volScalarField
        (
            IOobject
            (
                "XiEq",
                epsilon.time().timeName(),
                epsilon.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            epsilon.mesh(),
            dimensionedScalar(dimless, Zero)
        )
    );

    const objectRegistry& db = Su_.db();

    const volScalarField& b = db.lookupObject<volScalarField>("b");

    const volScalarField multiP1(0.0*pos(b - 0.99) + 1.0*neg(b - 0.99));
    const volScalarField multiP2(1.0*pos(b - 0.01) + 0.0*neg(b - 0.01));

    volScalarField& xieq = tXiEq.ref();

    forAll(xieq, celli)
    {
        scalar alpha;
        scalar beta;
        scalar gulderMa;

        if (Ma[celli]>= 0)
        {
            gulderMa =
                1.0
              + (0.1402 - 0.007*Ma[celli])
              * upBySu[celli]*sqrt(upBySu[celli]/K[celli]);

            regime[celli] = multiP1[celli]*multiP2[celli];

        }
        else
        {
            gulderMa =
                1.0
              + (0.005*Ma[celli]*Ma[celli]+0.01*Ma[celli] + 0.125)
              * upBySu[celli]*sqrt(upBySu[celli]/K[celli]);

             regime[celli] = 2*multiP1[celli]*multiP2[celli];
        }


        if (K[celli] < (lowK0_ + lowKg_*Ma[celli]) )
        {
            xieq[celli] = gulderMa;
        }
        else
        {
            if (Ma[celli] >= 0.0)
            {
                alpha = alphaCoefP_*(maLim_ - Ma[celli]);
                beta = betaCoefP_*(maLim_ - Ma[celli]);
                regime[celli] = 3*multiP1[celli]*multiP2[celli];
            }
            else
            {
                alpha = alphaCoefN_*(maLim1_ - Ma[celli]) ;
                beta = betaCoefN_*(maLim_ + Ma[celli]);
                regime[celli] = 4*multiP1[celli]*multiP2[celli];
            }
            xieq[celli] = XiEqCoef_*alpha*pow(K[celli], beta)*upBySu[celli];
        }

        if (Ma[celli] > -3.0 && Ma[celli] < 11.0)
        {
            scalar K0p8 = quenchCoef_*pow( Ma[celli] - quenchM_, quenchExp_);
            scalar quenchRate = quenchRate1_ + quenchRate2_*Ma[celli];
            if (K[celli] > (K0p8 - 0.223/quenchRate))
            {
                xieq[celli] *= 0.8*exp(-quenchRate*(K[celli] - K0p8));
                regime[celli] = 5*multiP1[celli]*multiP2[celli];
            }
        }
    }

    forAll(xieq.boundaryField(), patchi)
    {
        scalarField& xieqp = xieq.boundaryFieldRef()[patchi];
        const scalarField& Kp = K.boundaryField()[patchi];
        const scalarField& Map = Ma.boundaryField()[patchi];
        const scalarField& upBySup = upBySu.boundaryField()[patchi];

        forAll(xieqp, facei)
        {
            scalar alpha;
            scalar beta;

            if (Map[facei] > 0.0)
            {
                alpha = alphaCoefP_*(maLim_ - Map[facei]);
                beta = betaCoefP_*(maLim_ - Map[facei]);
            }
            else
            {
                alpha = alphaCoefN_*(maLim_ - Map[facei]);
                beta = betaCoefN_*(maLim_ + Map[facei]);
            }
            xieqp[facei] =
                XiEqCoef_*alpha*pow(Kp[facei], beta)*upBySup[facei];

            if (Map[facei] > -3.0 && Map[facei] < 11.0)
            {
                scalar K0p8 = quenchCoef_*pow(Map[facei] - quenchM_, quenchExp_);
                scalar quenchRate = quenchRate1_ + quenchRate2_*Ma[facei];

                if (Kp[facei] > (K0p8 - 0.223/quenchRate))
                {
                    xieqp[facei] *= 0.8*exp(-quenchRate*(Kp[facei] - K0p8));
                }
            }
            else
            {
                Info<<
                    "Markstein Number out of range for Quench Formulation" << endl;
            }
        }
    }

    return tXiEq;
}


bool Foam::XiEqModels::BLMgMaXiEq::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    return true;
}


// ************************************************************************* //
