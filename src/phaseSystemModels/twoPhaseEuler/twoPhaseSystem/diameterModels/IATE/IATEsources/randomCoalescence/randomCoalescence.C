/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "randomCoalescence.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(randomCoalescence, 0);
    addToRunTimeSelectionTable(IATEsource, randomCoalescence, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::randomCoalescence::
randomCoalescence
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Crc_("Crc", dimless, dict),
    C_("C", dimless, dict),
    alphaMax_("alphaMax", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::IATEsources::randomCoalescence::R() const
{
    auto tR = volScalarField::New
    (
        "R",
        IOobject::NO_REGISTER,
        iate_.phase().U().mesh(),
        dimensionedScalar(dimless/dimTime, Zero)
    );
    auto& R = tR.ref();

    scalar Crc = Crc_.value();
    scalar C = C_.value();
    scalar alphaMax = alphaMax_.value();
    volScalarField Ut(this->Ut());
    const volScalarField& alpha = phase();
    const volScalarField& kappai = iate_.kappai();
    scalar cbrtAlphaMax = cbrt(alphaMax);

    forAll(R, celli)
    {
        if (alpha[celli] < alphaMax - SMALL)
        {
            scalar cbrtAlphaMaxMAlpha = cbrtAlphaMax - cbrt(alpha[celli]);

            R[celli] =
                (-12)*phi()*kappai[celli]*alpha[celli]
               *Crc
               *Ut[celli]
               *(1 - exp(-C*cbrt(alpha[celli]*alphaMax)/cbrtAlphaMaxMAlpha))
               /(cbrtAlphaMax*cbrtAlphaMaxMAlpha);
        }
    }

    return tR;
}


// ************************************************************************* //
