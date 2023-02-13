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
    Foam::ChargeParcel

Group
    grpLagrangianIntermediateParcels

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes reacting multiphase parcel sub-models, plus:
    - electric charge

SourceFiles
    ChargeParcelI.H
    ChargeParcel.C
    ChargeParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ChargeParcel_H
#define Foam_ChargeParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType> class ChargeParcel;

template<class ParcelType>
Ostream& operator<<(Ostream&, const ChargeParcel<ParcelType>&);

/*---------------------------------------------------------------------------*\
                        Class ChargeParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ChargeParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    //- Class to hold charge particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);
    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


protected:

    // Protected data

        // Parcel properties

            //- Electric charge [A s] = [C]
            scalar eq_;


public:

    // Static data members

        //- Runtime type information
        TypeName("ChargeParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " eq"
        );

    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline ChargeParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell,
        //- searching for the rest of the required topology.
        //  Other properties are initialised as null.
        inline ChargeParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ChargeParcel
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
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const scalar eq0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ChargeParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        ChargeParcel(const ChargeParcel& p);

        //- Construct as a copy
        ChargeParcel
        (
            const ChargeParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ChargeParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ChargeParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //- parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ChargeParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<ChargeParcel<ParcelType>>
                (
                    new ChargeParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to electric charge
            inline scalar eq() const;


        // Edit

            //- Return access to electric charge
            inline scalar& eq();


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

            //- Read - no composition
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
            const ChargeParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChargeParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ChargeParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
