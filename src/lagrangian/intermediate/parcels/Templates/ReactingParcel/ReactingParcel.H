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
    Foam::ReactingParcel

Group
    grpLagrangianIntermediateParcels

Description
    Reacting parcel class with one/two-way coupling with the continuous
    phase.

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;



    //- Class to hold reacting parcel constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
    public:

            //- Type of activity coefficient models
            enum volumeUpdateType
            {
                mConstRho,
                mConstVol,
                mUpdateRhoAndVol,
                mUndefined
            };

    private:

        // Private Data

            //- Minimum pressure [Pa]
            demandDrivenEntry<scalar> pMin_;

            //- Method to update particle rho and diameter
            //- 0: constant rho
            //- 1: constant volume
            //- 2: recalculation of rho and volume based on the lost mass
            demandDrivenEntry<bool> constantVolume_;

            //- Method to update vol and rho
            demandDrivenEntry<label> volUpdateType_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the minimum pressure
            inline scalar pMin() const;

            //- Return const access to the constant volume flag
            inline bool constantVolume() const;

            //- Return const access to the constant volume flag
            inline label volUpdateType() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase pressure field
                autoPtr<interpolation<scalar>> pInterp_;


            // Cached continuous phase properties

                //- Pressure [Pa]
                scalar pc_;


    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions

            //- Return const access to the interpolator for continuous phase
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;

            //- Return the continuous phase pressure
            inline scalar pc() const;

            //- Access the continuous phase pressure
            inline scalar& pc();
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial mass [kg]
            scalar mass0_;

            //- Mass fractions of mixture []
            scalarField Y_;


    // Protected Member Functions

            //- Return change of volume due to mass exchange
            template<class TrackCloudType>
            scalar updatedDeltaVolume
            (
                TrackCloudType& cloud,
                const scalarField& dMass,
                const scalar p,
                const scalar T
            );


        //- Calculate Phase change
        template<class TrackCloudType>
        void calcPhaseChange
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const scalar rho,          // density
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YLiq,   // liquid component  mass fractions
            const scalarField& YSol,   // solid component mass fractions
            scalarField& dMassPC,      // mass transfer - local to parcel
            scalar& Sh,                // explicit parcel enthalpy source
            scalar& N,                 // flux of species emitted from parcel
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Update mass fraction
        scalar updateMassFraction
        (
            const scalar mass0,
            const scalarField& dMass,
            scalarField& Y
        ) const;


public:

    // Static data members

        //- Runtime type information
        TypeName("ReactingParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " mass0"
          + " nPhases(Y1..YN)"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ReactingParcel
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
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        ReactingParcel
        (
            const ReactingParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        ReactingParcel(const ReactingParcel& p);

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

            autoPtr<ReactingParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<ReactingParcel<ParcelType>>
                (
                    new ReactingParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial mass [kg]
            inline scalar mass0() const;

            //- Return const access to mass fractions of mixture []
            inline const scalarField& Y() const;

            //- Return const access to mass fractions of gases
            //  Note: for compatibilty only - returns Y()
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            //  Note: for compatibilty only - returns Y()
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            //  Note: for compatibilty only - returns Y()
            inline const scalarField& YSolid() const;


        // Edit

            //- Return access to initial mass [kg]
            inline scalar& mass0();

            //- Return access to mass fractions of mixture []
            inline scalarField& Y();


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

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
                scalar& Prs,
                scalar& kappas
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read - composition supplied
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write - composition supplied
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - no composition
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
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
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"
#include "ReactingParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
