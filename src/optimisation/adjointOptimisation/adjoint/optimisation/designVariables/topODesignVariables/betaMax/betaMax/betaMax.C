/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "betaMax.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(betaMax, 0);
    defineRunTimeSelectionTable(betaMax, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::betaMax::betaMax
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    value_(dict.getOrDefault<scalar>("betaMax", Zero))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::betaMax> Foam::betaMax::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.getOrDefault<word>("betaMaxType", "value"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    Info<< "betaMax type " << modelType << endl;

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "betaMaxType",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<betaMax>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::betaMax::value() const
{
    return value_;
}


// ************************************************************************* //
