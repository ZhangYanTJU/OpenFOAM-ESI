/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "waxSolventEvaporation.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "zeroField.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waxSolventEvaporation, 0);

addToRunTimeSelectionTable
(
    phaseChangeModel,
    waxSolventEvaporation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar waxSolventEvaporation::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waxSolventEvaporation::waxSolventEvaporation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    phaseChangeModel(typeName, film, dict),
    Wwax_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "Wwax"),
            film.regionMesh().time().constant(),
            film.regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        coeffDict_.get<scalar>("Wwax")
    ),
    Wsolvent_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "Wsolvent"),
            film.regionMesh().time().constant(),
            film.regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        coeffDict_.get<scalar>("Wsolvent")
    ),
    Ysolvent0_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "Ysolvent0"),
            film.regionMesh().time().constant(),
            film.regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    Ysolvent_
    (
        IOobject
        (
            IOobject::scopedName(typeName, "Ysolvent"),
            film.regionMesh().time().timeName(),
            film.regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        film.regionMesh()
    ),
    deltaMin_(coeffDict_.get<scalar>("deltaMin")),
    L_(coeffDict_.get<scalar>("L")),
    TbFactor_(coeffDict_.getOrDefault<scalar>("TbFactor", 1.1)),
    YInfZero_(coeffDict_.getOrDefault("YInfZero", false)),
    activityCoeff_
    (
        Function1<scalar>::New("activityCoeff", coeffDict_, &film.regionMesh())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

waxSolventEvaporation::~waxSolventEvaporation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class YInfType>
void waxSolventEvaporation::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy,
    const YInfType& YInf
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const volScalarField& delta = film.delta();
    const volScalarField& deltaRho = film.deltaRho();
    const surfaceScalarField& phi = film.phi();

    // Set local thermo properties
    const SLGThermo& thermo = film.thermo();
    const filmThermoModel& filmThermo = film.filmThermo();
    const label vapId = thermo.carrierId(filmThermo.name());

    // Retrieve fields from film model
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.T();
    const scalarField& hs = film.hs();
    const scalarField& rho = film.rho();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& magSf = film.magSf();
    const vectorField dU(film.UPrimary() - film.Us());
    const scalarField limMass
    (
        max(scalar(0), availableMass - deltaMin_*rho*magSf)
    );

    // Molecular weight of vapour [kg/kmol]
    const scalar Wvap = thermo.carrier().W(vapId);

    const scalar Wwax = Wwax_.value();
    const scalar Wsolvent = Wsolvent_.value();

    auto tevapRateCoeff = volScalarField::Internal::New
    (
        IOobject::scopedName(typeName, "evapRateCoeff"),
        film.regionMesh(),
        dimensionedScalar(dimDensity*dimVelocity, Zero)
    );
    auto& evapRateCoeff = tevapRateCoeff.ref();

    auto tevapRateInf = volScalarField::Internal::New
    (
        IOobject::scopedName(typeName, "evapRateInf"),
        film.regionMesh(),
        dimensionedScalar(dimDensity*dimVelocity, Zero)
    );
    auto& evapRateInf = tevapRateInf.ref();

    bool filmPresent = false;

    forAll(dMass, celli)
    {
        if (delta[celli] > deltaMin_)
        {
            filmPresent = true;

            const scalar Ysolvent = Ysolvent_[celli];

            // Molefraction of solvent in liquid film
            const scalar Xsolvent
            (
                Ysolvent*Wsolvent/((1 - Ysolvent)*Wwax + Ysolvent*Wsolvent)
            );

            // Primary region density [kg/m3]
            const scalar rhoInfc = rhoInf[celli];

            // Cell pressure [Pa]
            const scalar pc = pInf[celli];

            // Calculate the boiling temperature
            const scalar Tb = filmThermo.Tb(pc);

            // Local temperature - impose lower limit of 200 K for stability
            const scalar Tloc = min(TbFactor_*Tb, max(200.0, T[celli]));

            const scalar pPartialCoeff
            (
                filmThermo.pv(pc, Tloc)*activityCoeff_->value(Xsolvent)
            );

            scalar XsCoeff = pPartialCoeff/pc;

            // Vapour phase mole fraction of solvent at interface
            scalar Xs = XsCoeff*Xsolvent;

            if (Xs > 1)
            {
                WarningInFunction
                    << "Solvent vapour pressure > ambient pressure"
                    << endl;

                XsCoeff /= Xs;
                Xs = 1;
            }

            // Vapour phase mass fraction of solvent at the interface
            const scalar YsCoeff
            (
                XsCoeff/(XsCoeff*Xsolvent*Wsolvent + (1 - Xs)*Wvap)
            );

            // Primary region viscosity [Pa.s]
            const scalar muInfc = muInf[celli];

            // Reynolds number
            const scalar Re = rhoInfc*mag(dU[celli])*L_/muInfc;

            // Vapour diffusivity [m2/s]
            const scalar Dab = filmThermo.D(pc, Tloc);

            // Schmidt number
            const scalar Sc = muInfc/(rhoInfc*(Dab + ROOTVSMALL));

            // Sherwood number
            const scalar Sh = this->Sh(Re, Sc);

            // Mass transfer coefficient [m/s]
            evapRateCoeff[celli] = rhoInfc*Sh*Dab/(L_ + ROOTVSMALL);

            // Solvent mass transfer
            const scalar dm
            (
                max
                (
                    dt*magSf[celli]
                   *evapRateCoeff[celli]*(YsCoeff*Ysolvent - YInf[celli]),
                    0
                )
            );

            if (dm > limMass[celli])
            {
                evapRateCoeff[celli] *= limMass[celli]/dm;
            }

            evapRateInf[celli] = evapRateCoeff[celli]*YInf[celli];
            evapRateCoeff[celli] *= YsCoeff;

            // hVap[celli] = filmThermo.hl(pc, Tloc);
        }
    }

    const dimensionedScalar deltaRho0Bydt
    (
        "deltaRho0",
        deltaRho.dimensions()/dimTime,
        ROOTVSMALL/dt
    );

    volScalarField::Internal impingementRate
    (
        max
        (
           -film.rhoSp()(),
            dimensionedScalar(film.rhoSp().dimensions(), Zero)
        )
    );

    if (filmPresent)
    {
        // Solve for the solvent mass fraction
        fvScalarMatrix YsolventEqn
        (
            fvm::ddt(deltaRho, Ysolvent_)
          + fvm::div(phi, Ysolvent_)
         ==
            deltaRho0Bydt*Ysolvent_()

          + evapRateInf

            // Include the effect of the impinging droplets
            // added with Ysolvent = Ysolvent0
          + impingementRate*Ysolvent0_

          - fvm::Sp
            (
                deltaRho0Bydt
              + evapRateCoeff
              + film.rhoSp()()
              + impingementRate,
                Ysolvent_
            )
        );

        YsolventEqn.relax();
        YsolventEqn.solve();

        Ysolvent_.clamp_range(zero_one{});

        scalarField dm
        (
            dt*magSf*rhoInf*(evapRateCoeff*Ysolvent_ + evapRateInf)
        );

        dMass += dm;

        // Heat is assumed to be removed by heat-transfer to the wall
        // so the energy remains unchanged by the phase-change.
        dEnergy += dm*hs;

        // Latent heat [J/kg]
        // dEnergy += dm*(hs[celli] + hVap);
    }
}


void waxSolventEvaporation::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    if (YInfZero_)
    {
        correctModel(dt, availableMass, dMass, dEnergy, zeroField());
    }
    else
    {
        const thermoSingleLayer& film = filmType<thermoSingleLayer>();
        const label vapId = film.thermo().carrierId(film.filmThermo().name());
        const scalarField& YInf = film.YPrimary()[vapId];

        correctModel(dt, availableMass, dMass, dEnergy, YInf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
