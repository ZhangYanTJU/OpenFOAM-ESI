/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "parcelCloudModelList.H"

template<class ... Args>
void Foam::parcelCloudModelList::initialise(const Args& ... args)
{
    IOdictionary props
    (
        IOobject
        (
            "cloudProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ
        )
    );

    wordHashSet cloudTypes;
    if (props.readIfPresent("cloudTypes", cloudTypes))
    {
        setSize(cloudTypes.size());

        label i = 0;
        for (const auto& name : cloudTypes)
        {
            Info<< "creating cloud: " << name << endl;

            set(i++, parcelCloudModel::New(name, args ...));

            Info<< endl;
        }
    }
    else
    {
        setSize(1);
        set(0, parcelCloudModel::New("cloud", args ...));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parcelCloudModelList::parcelCloudModelList
(
    const dimensionedVector& g,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu
)
:
    PtrList<parcelCloudModel>(),
    mesh_(rho.mesh())
{
    initialise(g, rho, U, mu);
}


Foam::parcelCloudModelList::parcelCloudModelList
(
    const dimensionedVector& g,
    const volScalarField& rho,
    const volVectorField& U,
    const SLGThermo& slgThermo
)
:
    PtrList<parcelCloudModel>(),
    mesh_(rho.mesh())
{
    initialise(g, rho, U, slgThermo);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parcelCloudModelList::info()
{
    for (auto& c : *this)
    {
        c.info();
    }
}


void Foam::parcelCloudModelList::evolve()
{
    for (auto& c : *this)
    {
        c.evolveME();
    }
}


void Foam::parcelCloudModelList::storeGlobalPositions()
{
    for (auto& c : *this)
    {
        c.storeGlobalPositionsME();
    }
}


// ************************************************************************* //
