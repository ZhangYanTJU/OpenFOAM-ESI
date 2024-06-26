/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "GuldersEGR.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(GuldersEGR, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        GuldersEGR,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::GuldersEGR::GuldersEGR
(
    const dictionary& dict,
    const psiuReactionThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs").subDict(fuel_)),
    W_(coeffsDict_.get<scalar>("W")),
    eta_(coeffsDict_.get<scalar>("eta")),
    xi_(coeffsDict_.get<scalar>("xi")),
    f_(coeffsDict_.get<scalar>("f")),
    alpha_(coeffsDict_.get<scalar>("alpha")),
    beta_(coeffsDict_.get<scalar>("beta"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::GuldersEGR::~GuldersEGR()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::GuldersEGR::SuRef
(
    scalar phi
) const
{
    if (phi > SMALL)
    {
        return W_*pow(phi, eta_)*exp(-xi_*sqr(phi - 1.075));
    }
    else
    {
        return 0.0;
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    scalar p,
    scalar Tu,
    scalar phi,
    scalar Yres
) const
{
    static const scalar Tref = 300.0;
    static const scalar pRef = 1.013e5;

    return SuRef(phi)*pow((Tu/Tref), alpha_)*pow((p/pRef), beta_)*(1 - f_*Yres);
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    scalar phi
) const
{
    auto tSu0 = volScalarField::New
    (
        "Su0",
        IOobject::NO_REGISTER,
        p.mesh(),
        dimensionedScalar(dimVelocity, Zero)
    );
    auto& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi, 0.0);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        forAll(Su0Bf[patchi], facei)
        {
            Su0Bf[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    phi,
                    0.0
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& phi,
    const volScalarField& egr
) const
{
    auto tSu0 = volScalarField::New
    (
        "Su0",
        IOobject::NO_REGISTER,
        p.mesh(),
        dimensionedScalar(dimVelocity, Zero)
    );
    auto& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli], egr[celli]);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        forAll(Su0Bf[patchi], facei)
        {
            Su0Bf[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    phi.boundaryField()[patchi][facei],
                    egr.boundaryField()[patchi][facei]
                );
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::GuldersEGR::operator()() const
{
    if
    (
        psiuReactionThermo_.composition().contains("ft")
     && psiuReactionThermo_.composition().contains("egr")
    )
    {
        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            dimensionedScalar
            (
                "stoichiometricAirFuelMassRatio", dimless, psiuReactionThermo_
            )/
            (
                scalar(1)/psiuReactionThermo_.composition().Y("ft")
              - scalar(1)
            ),
            psiuReactionThermo_.composition().Y("egr")
        );
    }
    else
    {
        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            equivalenceRatio_
        );
    }
}


// ************************************************************************* //
