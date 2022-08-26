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
    Foam::ParcelCloudModel

Description
    Outermost template for parcel clouds. Adds the parcelCloudModel
    virtualisation layer and forwards the methods required by that layer.

SourceFiles
    ParcelCloudModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ParcelCloudModel_H
#define Foam_ParcelCloudModel_H

#include "parcelCloudModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
class ParcelCloudModel
:
    public parcelCloudModel,
    public CloudType
{
public:

    // Constructors

        //- Inherit base constructors
        using CloudType::CloudType;


    // Member Functions

        //- Evolve the cloud
        virtual void evolveME() final
        {
            // CloudType::evolve();
            if (this->solution().canEvolve())
            {
                typedef typename CloudType::particleType parcelType;

                typename parcelType::trackingData td(*this);

                this->solve(*this, td);
            }
        }

        //- Call this before a topology change. Stores the particles global
        //- positions in the database for use during mapping.
        virtual void storeGlobalPositionsME() const final
        {
            CloudType::storeGlobalPositions();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
