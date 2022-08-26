/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "parcelCloudModel.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::parcelCloudModel> Foam::parcelCloudModel::New
(
    const word& cloudName,
    const dimensionedVector& g,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu
)
{
    const fvMesh& mesh = rho.mesh();

    IOdictionary dict
    (
        IOobject
        (
            cloudName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType(dict.get<word>("type"));

    Info<< "Selecting cloud type " << modelType << endl;

    auto* ctorPtr = componentsConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *componentsConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<parcelCloudModel>(ctorPtr(cloudName, g, rho, U, mu));
}


Foam::autoPtr<Foam::parcelCloudModel> Foam::parcelCloudModel::New
(
    const word& cloudName,
    const dimensionedVector& g,
    const volScalarField& rho,
    const volVectorField& U,
    const SLGThermo& thermo
)
{
    const fvMesh& mesh = rho.mesh();

    IOdictionary dict
    (
        IOobject
        (
            cloudName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType(dict.get<word>("type"));

    Info<< "Selecting cloud type " << modelType << endl;

    auto* ctorPtr = thermoConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *thermoConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<parcelCloudModel>(ctorPtr(cloudName, g, rho, U, thermo));
}


// ************************************************************************* //