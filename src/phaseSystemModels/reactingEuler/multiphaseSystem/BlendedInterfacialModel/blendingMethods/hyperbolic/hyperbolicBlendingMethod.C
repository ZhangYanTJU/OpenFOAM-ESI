/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "hyperbolicBlendingMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(hyperbolic, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        hyperbolic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::hyperbolic::hyperbolic
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    transitionAlphaScale_("transitionAlphaScale", dimless, dict)
{
    for (const word& phaseName : phaseNames)
    {
        minContinuousAlpha_.insert
        (
            phaseName,
            dimensionedScalar
            (
                IOobject::groupName("minContinuousAlpha", phaseName),
                dimless,
                dict
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase2 - minContinuousAlpha_[phase2.name()])
            )
        )/2;
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase1 - minContinuousAlpha_[phase1.name()])
            )
        )/2;
}


// ************************************************************************* //
