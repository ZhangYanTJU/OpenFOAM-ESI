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
    Foam::molecule

Description
    Foam::molecule

SourceFiles
    moleculeI.H
    molecule.C
    moleculeIO.C

\*---------------------------------------------------------------------------*/

#ifndef molecule_H
#define molecule_H

#include "particle.H"
#include "IOstream.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class moleculeCloud;


// Forward declaration of friend functions and operators

class molecule;

Ostream& operator<<(Ostream&, const molecule&);


/*---------------------------------------------------------------------------*\
                           Class molecule Declaration
\*---------------------------------------------------------------------------*/

class molecule
:
    public particle
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


    // Values of special that are less than zero are for built-in functionality.
    // Values greater than zero are user specifiable/expandable
    // (i.e. test special_ >= SPECIAL_USER)

    enum specialTypes
    {
        SPECIAL_TETHERED = -1,
        SPECIAL_FROZEN   = -2,
        NOT_SPECIAL      = 0,
        SPECIAL_USER     = 1
    };


    //- Class to hold molecule constant properties
    class constantProperties
    {
        // Private data

            Field<vector> siteReferencePositions_;

            List<scalar> siteMasses_;

            List<scalar> siteCharges_;

            List<label> siteIds_;

            List<bool> pairPotentialSites_;

            List<bool> electrostaticSites_;

            diagTensor momentOfInertia_;

            scalar mass_;


        // Private Member Functions

            void checkSiteListSizes() const;

            void setInteracionSiteBools
            (
                const List<word>& siteIds,
                const List<word>& pairPotSiteIds
            );

            bool linearMoleculeTest() const;


    public:

        inline constantProperties();

        //- Construct from dictionary
        inline constantProperties(const dictionary& dict);

        // Member functions

        inline const Field<vector>& siteReferencePositions() const;

        inline const List<scalar>& siteMasses() const;

        inline const List<scalar>& siteCharges() const;

        inline const List<label>& siteIds() const;

        inline List<label>& siteIds();

        inline const List<bool>& pairPotentialSites() const;

        inline bool pairPotentialSite(label sId) const;

        inline const List<bool>& electrostaticSites() const;

        inline bool electrostaticSite(label sId) const;

        inline const diagTensor& momentOfInertia() const;

        inline bool linearMolecule() const;

        inline bool pointMolecule() const;

        inline label degreesOfFreedom() const;

        inline scalar mass() const;

        inline label nSites() const;
    };


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // label specifying which part of the integration algorithm is taking
        label part_;


    public:

        // Constructors

            trackingData(moleculeCloud& cloud, label part)
            :
                particle::trackingData(cloud),
                part_(part)
            {}

        // Member functions

            inline label part() const
            {
                return part_;
            }
    };


private:

    // Private data

        tensor Q_;

        vector v_;

        vector a_;

        vector pi_;

        vector tau_;

        vector specialPosition_;

        scalar potentialEnergy_;

        // - r_ij f_ij, stress dyad
        tensor rf_;

        label special_;

        label id_;

        List<vector> siteForces_;

        List<vector> sitePositions_;


    // Private Member Functions

        tensor rotationTensorX(scalar deltaT) const;

        tensor rotationTensorY(scalar deltaT) const;

        tensor rotationTensorZ(scalar deltaT) const;


public:

    friend class Cloud<molecule>;

    // Constructors

        //- Construct from components
        inline molecule
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
            const constantProperties& constProps,
            const label special,
            const label id
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        inline molecule
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
            const constantProperties& constProps,
            const label special,
            const label id
        );

        //- Construct from Istream
        molecule
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Return a clone
        virtual autoPtr<particle> clone() const
        {
            return particle::Clone(*this);
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

            autoPtr<molecule> operator()(Istream& is) const
            {
                return autoPtr<molecule>::New(mesh_, is, true);
            }
        };


    // Member Functions

        // Tracking

            bool move(moleculeCloud&, trackingData&, const scalar trackTime);

            virtual void transformProperties(const tensor& T);

            virtual void transformProperties(const vector& separation);

            void setSitePositions(const constantProperties& constProps);

            void setSiteSizes(label size);


        // Access

            inline const tensor& Q() const;
            inline tensor& Q();

            inline const vector& v() const;
            inline vector& v();

            inline const vector& a() const;
            inline vector& a();

            inline const vector& pi() const;
            inline vector& pi();

            inline const vector& tau() const;
            inline vector& tau();

            inline const List<vector>& siteForces() const;
            inline List<vector>& siteForces();

            inline const List<vector>& sitePositions() const;
            inline List<vector>& sitePositions();

            inline const vector& specialPosition() const;
            inline vector& specialPosition();

            inline scalar potentialEnergy() const;
            inline scalar& potentialEnergy();

            inline const tensor& rf() const;
            inline tensor& rf();

            inline label special() const;

            inline bool tethered() const;

            inline label id() const;


    // Member Operators

        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
        bool hitPatch(moleculeCloud& cloud, trackingData& td);

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch(moleculeCloud& cloud, trackingData& td);

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch(moleculeCloud& cloud, trackingData& td);


    // I-O

        static void readFields(Cloud<molecule>& mC);

        static void writeFields(const Cloud<molecule>& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const molecule&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
