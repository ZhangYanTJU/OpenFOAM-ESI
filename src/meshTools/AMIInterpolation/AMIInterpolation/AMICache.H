/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::AMICache

Description
    Provides caching of weights and addressing to AMIInterpolation

SourceFiles
    AMICache.C

SeeAlso
    Foam::AMIInterpolation.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_AMICache_H
#define Foam_AMICache_H

#include "cylindricalCS.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class AMIInterpolation;

/*---------------------------------------------------------------------------*\
                          Class AMICache Declaration
\*---------------------------------------------------------------------------*/

class AMICache
{
public:

    // Public Data Types

        //- Tolerance used when caching the AMI to identify e.g. if the
        //- current rotation angle has already been captured
        static scalar cacheThetaTolerance_;

        static int debug;

private:

    // Private Data

        //- Cache size
        label size_;

        //- Axis of rotation for rotational cyclics
        vector rotationAxis_;

        //- Point on axis of rotation for rotational cyclics
        point rotationCentre_;

        //- Flag to indicate that the cache is complete
        bool complete_;

        //- Cache index 0 in lists
        label index0_;

        //- Cache index 1 in lists
        label index1_;

        //- Interpolation weight between cache indices 0 and 1
        scalar interpWeight_;

        //- Local co-ordinate system
        autoPtr<coordSystem::cylindrical> coordSysPtr_;

        //- List of cached angle snapshots
        List<scalar> theta_;

        //- List of source addresses
        List<labelListList> cachedSrcAddress_;

        //- List of source weights
        List<scalarListList> cachedSrcWeights_;

        //- List of source weights sums
        List<scalarField> cachedSrcWeightsSum_;

        //- List of source parallel maps
        List<autoPtr<mapDistribute>> cachedSrcMapPtr_;

        //- List of target addresses
        List<labelListList> cachedTgtAddress_;

        //- List of target weights
        List<scalarListList> cachedTgtWeights_;

        //- List of target weights sums
        List<scalarField> cachedTgtWeightsSum_;

        //- List of target parallel maps
        List<autoPtr<mapDistribute>> cachedTgtMapPtr_;


        //- Flag to indicate interpolation direction
        mutable bool toSource_;


    // Private Member Functions

        //- Get rotation angle for point using local co-ordinate system
        scalar getRotationAngle(const point& globalPoint) const;


public:

    // Constructors

        //- Null constructor
        AMICache(const bool toSource = true);

        //- Construct from dictionary
        AMICache(const dictionary& dict, const bool toSource = true);

        //- Construct as copy
        AMICache(const AMICache& cache);

        //- Construct from agglomeration of AMIInterpolation. Agglomeration
        //- passed in as new coarse size and addressing from fine from coarse
        AMICache
        (
            const AMICache& cache,
            const AMIInterpolation& fineAMI,
            const labelList& sourceRestrictAddressing,
            const labelList& targetRestrictAddressing
        );

        //- Construct from stream
        AMICache(Istream& is);


    // Member Functions

        //- Return true if cache is active
        constexpr bool active() const noexcept { return size_ > 0; }

        //- Return cache size
        constexpr label size() const noexcept { return size_; }

        //- Return true if cache is complete
        constexpr bool complete() const noexcept { return complete_; }

        //- Return cache lower bound index
        constexpr label index0() const noexcept { return index0_; }

        //- Return cache upper bound index
        constexpr label index1() const noexcept { return index1_; }

        //- Return cache interpolation weight
        constexpr label weight() const noexcept { return interpWeight_; }

        //- Return list of cached rotation angles
        const List<scalar>& theta() const noexcept { return theta_; }

        //- Return List of source addresses
        const List<labelListList>& cachedSrcAddress() const noexcept
        {
            return cachedSrcAddress_;
        }

        //- Return List of source weights
        const List<scalarListList>& cachedSrcWeights() const noexcept
        {
            return cachedSrcWeights_;
        }

        //- Return List of source weights sums
        const List<scalarField>& cachedSrcWeightsSum() const noexcept
        {
            return cachedSrcWeightsSum_;
        }

        //- Return List of source parallel maps
        const List<autoPtr<mapDistribute>>& cachedSrcMapPtr() const noexcept
        {
            return cachedSrcMapPtr_;
        }

        //- Return List of target addresses
        const List<labelListList>& cachedTgtAddress() const noexcept
        {
            return cachedTgtAddress_;
        }

        //- Return List of target weights
        const List<scalarListList>& cachedTgtWeights() const noexcept
        {
            return cachedTgtWeights_;
        }

        //- Return List of target weights sums
        const List<scalarField>& cachedTgtWeightsSum() const noexcept
        {
            return cachedTgtWeightsSum_;
        }

        //- Return List of target parallel maps
        const List<autoPtr<mapDistribute>>& cachedTgtMapPtr() const noexcept
        {
            return cachedTgtMapPtr_;
        }

        //- Apply cached evaluation based on user supplied evaluation function
        template<class Type, class EvalFunction>
        bool apply(List<Type>& result, const EvalFunction& eval) const;

        //- Flag that lower bound is applicable
        constexpr bool applyLower() const noexcept
        {
            return index0_ != -1 && index1_ == -1;
        }

        //- Flag that upper bound is applicable
        constexpr bool applyUpper() const noexcept
        {
            return index0_ == -1 && index1_ != -1;
        }

        //- Flag that interpolation is applicable
        constexpr bool applyInterpolate() const noexcept
        {
            return index0_ != -1 && index1_ != -1;
        }

        //- Set the interpolation direction
        constexpr bool setDirection(bool toSource) const noexcept
        {
            toSource_ = toSource;
            return toSource_;
        }

        //- Restore AMI weights and addressing from the cache
        bool restoreCache(const point& globalPoint);

        //- Add AMI weights and addressing to the cache
        void addToCache
        (
            const AMIInterpolation& ami,
            const point& globalPoint
        );

        //- Check cache index is valid
        void checkIndex(const label index) const
        {
            if ((index < 0) || (index > size_-1))
            {
                FatalErrorInFunction
                    << "Supplied out of bounds: " << index << "/"
                    << size_ << abort(FatalError);
            }
        }


        // Helper functions to retrieve cached values at index0 and index1
        // Note: uses interpolation direction toSource_
        #undef  defineMethods01
        #define defineMethods01(Src, Tgt, idx)                                 \
            const labelListList& c##Src##Address##idx() const                  \
            {                                                                  \
                checkIndex(index##idx##_);                                     \
                return toSource_ ?                                             \
                    cached##Src##Address_[index##idx##_]                       \
                  : cached##Tgt##Address_[index##idx##_];                      \
            }                                                                  \
            const scalarListList& c##Src##Weights##idx() const                 \
            {                                                                  \
                checkIndex(index##idx##_);                                     \
                return toSource_ ?                                             \
                    cached##Src##Weights_[index##idx##_]                       \
                  : cached##Tgt##Weights_[index##idx##_];                      \
            }                                                                  \
            const scalarField& c##Src##WeightsSum##idx() const                 \
            {                                                                  \
                checkIndex(index##idx##_);                                     \
                return toSource_ ?                                             \
                    cached##Src##WeightsSum_[index##idx##_]                    \
                  : cached##Tgt##WeightsSum_[index##idx##_];                   \
            }                                                                  \
            const autoPtr<mapDistribute>& c##Src##MapPtr##idx() const          \
            {                                                                  \
                checkIndex(index##idx##_);                                     \
                return toSource_ ?                                             \
                    cached##Src##MapPtr_[index##idx##_]                        \
                  : cached##Tgt##MapPtr_[index##idx##_];                       \
            }

        defineMethods01(Src, Tgt, 0)
        defineMethods01(Src, Tgt, 1)
        defineMethods01(Tgt, Src, 0)
        defineMethods01(Tgt, Src, 1)


        // Helper functions to retrieve cached values at supplied index
        // Note: uses interpolation direction toSource_
        #undef  defineMethodsIndex
        #define defineMethodsIndex(Src, Tgt)                                   \
            const labelListList& c##Src##Address(const label index) const      \
            {                                                                  \
                checkIndex(index);                                             \
                return toSource_ ?                                             \
                    cached##Src##Address_[index]                               \
                  : cached##Tgt##Address_[index];                              \
            }                                                                  \
            const scalarListList& c##Src##Weights(const label index) const     \
            {                                                                  \
                checkIndex(index);                                             \
                return toSource_ ?                                             \
                    cached##Src##Weights_[index]                               \
                  : cached##Tgt##Weights_[index];                              \
            }                                                                  \
            const scalarField& c##Src##WeightsSum(const label index) const     \
            {                                                                  \
                checkIndex(index);                                             \
                return toSource_ ?                                             \
                    cached##Src##WeightsSum_[index]                            \
                  : cached##Tgt##WeightsSum_[index];                           \
            }                                                                  \
            const autoPtr<mapDistribute>& c##Src##MapPtr(const label index)    \
            const                                                              \
            {                                                                  \
                checkIndex(index);                                             \
                return toSource_ ?                                             \
                    cached##Src##MapPtr_[index]                                \
                  : cached##Tgt##MapPtr_[index];                               \
            }

        defineMethodsIndex(Src, Tgt)
        defineMethodsIndex(Tgt, Src)


    // I-O

        //- Write AMI as a dictionary
        void write(Ostream& os) const;

        //- Write AMI raw
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMICacheTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
