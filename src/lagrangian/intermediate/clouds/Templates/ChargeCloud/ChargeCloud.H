/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::ChargeCloud

Group
    grpLagrangianIntermediateClouds

Description
    Templated base class for electric-charge cloud

    - Adds to multiphase reacting cloud
      - electric charge

SourceFiles
    ChargeCloudI.H
    ChargeCloud.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ChargeCloud_H
#define Foam_ChargeCloud_H

#include "className.H"
#include "autoPtr.H"
#include "Cloud.H"
#include "volFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class SLGThermo;

TemplateName(ChargeCloud);

/*---------------------------------------------------------------------------*\
                        Class ChargeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ChargeCloud
:
    public CloudType,
    public ChargeCloudName
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ChargeCloud<CloudType> chargeCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<ChargeCloud<CloudType>> cloudCopyPtr_;


    // Private Member Functions

        //- No copy construct
        ChargeCloud(const ChargeCloud&) = delete;

        //- No copy assignment
        void operator=(const ChargeCloud&) = delete;


protected:

    // Protected Data

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud Evolution Functions

            //- Reset state of cloud
            void cloudReset(ChargeCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        ChargeCloud
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        ChargeCloud
        (
            ChargeCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        ChargeCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ChargeCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ChargeCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ChargeCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ChargeCloud() = default;


    // Member Functions

        // Access

        //- Return a reference to the cloud copy
        inline const ChargeCloud& cloudCopy() const;

        //- Return the constant properties
        inline const typename parcelType::constantProperties&
            constProps() const;

        //- Return access to the constant properties
        inline typename parcelType::constantProperties& constProps();


        // Sub-models


        // Cloud Evolution

        //- Set parcel thermo properties
        void setParcelThermoProperties
        (
            parcelType& parcel,
            const scalar lagrangianDt
        );

        //- Check parcel properties
        void checkParcelProperties
        (
            parcelType& parcel,
            const scalar lagrangianDt,
            const bool fullyDescribed
        );

        //- Store the current cloud state
        void storeState();

        //- Reset the current cloud to the previously stored state
        void restoreState();

        //- Reset the cloud source terms
        void resetSourceTerms();

        //- Evolve the cloud
        void evolve();


        // Mapping

        //- Remap the cells of particles corresponding to the
        //- mesh topology change with a default tracking data object
        virtual void autoMap(const mapPolyMesh&);


        // I-O

        //- Print cloud information
        void info();

        //- Write the field data for the cloud
        // virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChargeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChargeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
