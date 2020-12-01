/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::RecycleInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Parcel recycle patch interaction model.

    Parcels that hit an 'outflow' patch are recycled to an 'inflow' patch,
    with optional:
    - recycle fraction < 1

    Example usage:
    \verbatim
    RecycleInteractionCoeffs
    {
        recyclePatches
        (
            (outlet1 inlet1)
            (outlet2 inlet2)
        );

        recycleFraction     0.8; // [0-1]
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef RecycleInteraction_H
#define RecycleInteraction_H

#include "patchInjectionBase.H"
#include "IDLList.H"
#include "PtrList.H"
#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class RecycleInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RecycleInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    //typedef Tuple2<scalar, typename CloudType::parcelType> storedParcel;

    typedef typename CloudType::parcelType parcelType;

    // Protected data

        // Reference to mesh
        const fvMesh& mesh_;

        List<Tuple2<word, word>> recyclePatches_;

        List<Tuple2<label, label>> recyclePatchesIds_;

        List<IDLList<parcelType>> recycledParcels_;

        const scalar recycleFraction_;


        // Bookkeeping for particle fates

            //- Number of parcels removed
            List<List<label>> nRemoved_;

            //- Mass of parcels removed
            List<List<scalar>> massRemoved_;

            //- Number of parcels injected
            List<List<label>> nInjected_;

            //- Mass of parcels injected
            List<List<scalar>> massInjected_;


        //- Flag to output escaped/mass particles sorted by injectorID
        bool outputByInjectorId_;

        //- Injector ID to local index map
        Map<label> injIdToIndex_;

        PtrList<patchInjectionBase> injectionPatchPtr_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("recycleInteraction");


    // Constructors

        //- Construct from dictionary
        RecycleInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        RecycleInteraction(const RecycleInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new RecycleInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RecycleInteraction() = default;


    // Member Functions

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );

        //- Post-evolve hook
        virtual void postEvolve();

        //- Write patch interaction info to stream
        virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RecycleInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
