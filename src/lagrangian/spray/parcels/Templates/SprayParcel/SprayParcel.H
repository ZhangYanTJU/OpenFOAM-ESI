/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::SprayParcel

Description
    Reacting spray parcel, with added functionality for atomization and breakup

\*---------------------------------------------------------------------------*/

#ifndef SprayParcel_H
#define SprayParcel_H

#include "particle.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class SprayParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const SprayParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                          Class SprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class SprayParcel
:
    public ParcelType
{
public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Particle initial surface tension [N/m]
            demandDrivenEntry<scalar> sigma0_;

            //- Particle initial dynamic viscosity [Pa.s]
            demandDrivenEntry<scalar> mu0_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);

            //- Construct from components
            constantProperties
            (
                const label parcelTypeId,
                const scalar rhoMin,
                const scalar rho0,
                const scalar minParcelMass,
                const scalar youngsModulus,
                const scalar poissonsRatio,
                const scalar T0,
                const scalar TMin,
                const scalar TMax,
                const scalar Cp0,
                const scalar epsilon0,
                const scalar f0,
                const scalar Pr,
                const scalar pMin,
                const bool constantVolume,
                const scalar sigma0,
                const scalar mu0
            );


        // Access

            //- Return const access to the initial surface tension
            inline scalar sigma0() const;

            //- Return const access to the initial dynamic viscosity
            inline scalar mu0() const;
    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


protected:

    // Protected data

        // Spray parcel properties

            //- Initial droplet diameter
            scalar d0_;

            //- Injection position
            vector position0_;

            //- Liquid surface tension [N/m]
            scalar sigma_;

            //- Liquid dynamic viscosity [Pa.s]
            scalar mu_;

            //- Part of liquid core ( >0.5=liquid, <0.5=droplet )
            scalar liquidCore_;

            //- Index for KH Breakup
            scalar KHindex_;

            //- Spherical deviation
            scalar y_;

            //- Rate of change of spherical deviation
            scalar yDot_;

            //- Characteristic time (used in atomization and/or breakup model)
            scalar tc_;

            //- Stripped parcel mass due to breakup
            scalar ms_;

            //- Injected from injector (needed e.g. for calculating distance
            //  from injector)
            scalar injector_;

            //- Momentum relaxation time (needed for calculating parcel acc.)
            scalar tMom_;

            //- Passive scalar (extra variable to be defined by user)
            scalar user_;


public:

    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;

        //- Runtime type information
        TypeName("SprayParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " d0"
          + " position0"
          + " sigma"
          + " mu"
          + " liquidCore"
          + " KHindex"
          + " y"
          + " yDot"
          + " tc"
          + " ms"
          + " injector"
          + " tMom"
          + " user"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline SprayParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline SprayParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline SprayParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalarField& Y0,
            const scalar liquidCore,
            const scalar KHindex,
            const scalar y,
            const scalar yDot,
            const scalar tc,
            const scalar ms,
            const scalar injector,
            const scalar tMom,
            const scalar user,
            const typename ParcelType::constantProperties& constProps
        );

        //- Construct from Istream
        SprayParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        SprayParcel
        (
            const SprayParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        SprayParcel(const SprayParcel& p);

        //- Return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return particle::Clone(*this);
        }

        //- Return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return particle::Clone(*this, mesh);
        }

        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<SprayParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<SprayParcel<ParcelType>>
                (
                    new SprayParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial droplet diameter
            inline scalar d0() const;

            //- Return const access to initial droplet position
            inline const vector& position0() const;

            //- Return const access to the liquid surface tension
            inline scalar sigma() const;

            //- Return const access to the liquid dynamic viscosity
            inline scalar mu() const;

            //- Return const access to liquid core
            inline scalar liquidCore() const;

            //- Return const access to Kelvin-Helmholtz breakup index
            inline scalar KHindex() const;

            //- Return const access to spherical deviation
            inline scalar y() const;

            //- Return const access to rate of change of spherical deviation
            inline scalar yDot() const;

            //- Return const access to atomization characteristic time
            inline scalar tc() const;

            //- Return const access to stripped parcel mass
            inline scalar ms() const;

            //- Return const access to injector id
            inline scalar injector() const;

            //- Return const access to momentum relaxation time
            inline scalar tMom() const;

            //- Return const access to passive user scalar
            inline scalar user() const;


        // Edit

            //- Return access to initial droplet diameter
            inline scalar& d0();

            //- Return access to initial droplet position
            inline vector& position0();

            //- Return access to the liquid surface tension
            inline scalar& sigma();

            //- Return access to the liquid dynamic viscosity
            inline scalar& mu();

            //- Return access to liquid core
            inline scalar& liquidCore();

            //- Return access to Kelvin-Helmholtz breakup index
            inline scalar& KHindex();

            //- Return access to spherical deviation
            inline scalar& y();

            //- Return access to rate of change of spherical deviation
            inline scalar& yDot();

            //- Return access to atomization characteristic time
            inline scalar& tc();

            //- Return access to stripped parcel mass
            inline scalar& ms();

            //- Return access to injector id
            inline scalar& injector();

            //- Return access to momentum relaxation time
            inline scalar& tMom();

            //- Return access to passive user scalar
            inline scalar& user();


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct parcel properties according to atomization model
            template<class TrackCloudType>
            void calcAtomization
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct parcel properties according to breakup model
            template<class TrackCloudType>
            void calcBreakup
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct surface values due to emitted species
            template<class TrackCloudType>
            void correctSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappa
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Calculate the chi-factor for flash-boiling for the
            //  atomization model
            template<class TrackCloudType>
            scalar chi
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalarField& X
            ) const;

            //- Solve the TAB equation
            template<class TrackCloudType>
            void solveTABEq
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly
            ) const;

            //- Read particle fields as objects from the obr registry
            //  - no composition
            template<class CloudType>
            static void readObjects
            (
                CloudType& c,
                const objectRegistry& obr
            );

            //- Read particle fields as objects from the obr registry
            template<class CloudType, class CompositionType>
            static void readObjects
            (
                CloudType& c,
                const CompositionType& compModel,
                const objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            //  - no composition
            template<class CloudType>
            static void writeObjects
            (
                const CloudType& c,
                objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            template<class CloudType, class CompositionType>
            static void writeObjects
            (
                const CloudType& c,
                const CompositionType& compModel,
                objectRegistry& obr
            );


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const SprayParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SprayParcel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
