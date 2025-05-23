/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::MPPICParcel

Group
    grpLagrangianIntermediateParcels

Description
    Wrapper around kinematic parcel types to add MPPIC modelling

SourceFiles
    MPPICParcelI.H
    MPPICParcelTrackingDataI.H
    MPPICParcel.C
    MPPICParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef MPPICParcel_H
#define MPPICParcel_H

#include "particle.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class ParcelType>
class MPPICParcel;

template<class Type>
class AveragingMethod;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const MPPICParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class MPPICParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class MPPICParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    class trackingData
    :
        public ParcelType::trackingData
    {

    public:

        enum trackPart
        {
            tpLinearTrack,
            tpDampingNoTrack,
            tpPackingNoTrack,
            tpCorrectTrack,
        };


    private:

        // Private data

            // MPPIC Averages

                //- Volume average
                autoPtr<AveragingMethod<scalar>> volumeAverage_;

                //- Radius average [ volume^(1/3) ]
                autoPtr<AveragingMethod<scalar>> radiusAverage_;

                //- Density average
                autoPtr<AveragingMethod<scalar>> rhoAverage_;

                //- Velocity average
                autoPtr<AveragingMethod<vector>> uAverage_;

                //- Magnitude velocity squared average
                autoPtr<AveragingMethod<scalar>> uSqrAverage_;

                //- Frequency average
                autoPtr<AveragingMethod<scalar>> frequencyAverage_;

                //- Mass average
                autoPtr<AveragingMethod<scalar>> massAverage_;


            //- Label specifying the current part of the tracking process
            trackPart part_;


    public:

        //- Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = tpLinearTrack
            );


        //- Update the MPPIC averages
        template<class TrackCloudType>
        inline void updateAverages(const TrackCloudType& cloud);


        //- Access

            //- Const access to the tracking part label
            inline trackPart part() const;

            //- Non const access to the tracking part label
            inline trackPart& part();
    };


protected:

    // Protected data

        //- Velocity correction due to collisions [m/s]
        vector UCorrect_;


public:

    // Static data members

        //- Runtime type information
        TypeName("MPPICParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " (UCorrectx UCorrecty UCorrectz)"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline MPPICParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline MPPICParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline MPPICParcel
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
            const vector& UCorrect0,
            const typename ParcelType::constantProperties& constProps
        );

        //- Construct from Istream
        MPPICParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        MPPICParcel(const MPPICParcel& p);

        //- Construct as a copy
        MPPICParcel(const MPPICParcel& p, const polyMesh& mesh);

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

            autoPtr<MPPICParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<MPPICParcel<ParcelType>>
                (
                    new MPPICParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to correction velocity
            inline const vector& UCorrect() const;

            //- Return access to correction velocity
            inline vector& UCorrect();


        // Tracking

            //- Move the parcel
            template<class TrackCloudType>
            bool move
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar trackTime
            );


    // Friend Functions

        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
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
            template<class CloudType>
            static void readObjects(CloudType& c, const objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            template<class CloudType>
            static void writeObjects(const CloudType& c, objectRegistry& obr);


        // Ostream operator

            friend Ostream& operator<< <ParcelType>
            (
                Ostream&,
                const MPPICParcel<ParcelType>&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MPPICParcelI.H"
#include "MPPICParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MPPICParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
