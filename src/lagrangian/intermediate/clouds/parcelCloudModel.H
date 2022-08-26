/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenFOAM Foundation
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

Class
    Foam::parcelCloudModel

Description
    Virtual abstract base class for parcel clouds. As baseCloudInterface but
    with additional virtualisation of the evolve method, plus some additional
    methods that are defined below the parcel-cloud layer (i.e., in Cloud).
    These methods are implemented by forwarding in the ParcelCloudModel class.

    This class contains the selection table and represents the high level
    interface used by a solver.

SourceFiles
    parcelCloudModel.C
    parcelCloudModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parcelCloudModel_H
#define Foam_parcelCloudModel_H

#include "dimensionedVector.H"
#include "baseCloudInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SLGThermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class parcelCloudModel
:
    virtual public baseCloudInterface
{
public:

    //- Runtime type information
    TypeName("parcelCloudModel");

    //- Null constructor
    parcelCloudModel()
    {}

    //- Destructor
    virtual ~parcelCloudModel() = default;


    declareRunTimeSelectionTable
    (
        autoPtr,
        parcelCloudModel,
        components,
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu
        ),
        (cloudName, g, rho, U, mu)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        parcelCloudModel,
        thermo,
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const SLGThermo& thermo
        ),
        (cloudName, g, rho, U, thermo)
    );

    // Selectors

        static autoPtr<parcelCloudModel> New
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu
        );

        static autoPtr<parcelCloudModel> New
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const SLGThermo& thermo
        );


    // Member Functions

        virtual void storeGlobalPositionsME() const = 0;

        virtual void evolveME() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
