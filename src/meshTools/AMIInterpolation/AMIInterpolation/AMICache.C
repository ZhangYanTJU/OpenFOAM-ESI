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

\*---------------------------------------------------------------------------*/

#include "AMICache.H"
#include "AMIInterpolation.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::AMICache::cacheThetaTolerance_ = 1e-8;

int Foam::AMICache::debug = 0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::AMICache::getRotationAngle(const point& globalPoint) const
{
    if (!coordSysPtr_)
    {
        FatalErrorInFunction
            << "No co-ordinate system available for theta evaluation"
            << abort(FatalError);
    }


    scalar theta = coordSysPtr_->localPosition(globalPoint).y();

    // Ensure 0 < theta < 2pi
    if (mag(theta) < cacheThetaTolerance_)
    {
        theta = 0;
    }
    else if (theta < 0)
    {
        theta += constant::mathematical::twoPi;
    }

    return theta;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AMICache::AMICache(const dictionary& dict, const bool toSource)
:
    size_(dict.getOrDefault<label>("cacheSize", 0)),
    complete_(false),
    index0_(-1),
    index1_(-1),
    interpWeight_(0),
    coordSysPtr_(nullptr),
    theta_(),
    cachedSrcAddress_(),
    cachedSrcWeights_(),
    cachedSrcWeightsSum_(),
    cachedSrcMapPtr_(),
    cachedTgtAddress_(),
    cachedTgtWeights_(),
    cachedTgtWeightsSum_(),
    cachedTgtMapPtr_(),
    toSource_(toSource)
{
    if (size_ != 0)
    {
        theta_.resize(size_, GREAT);
        cachedSrcAddress_.resize(size_);
        cachedSrcWeights_.resize(size_);
        cachedSrcWeightsSum_.resize(size_);
        cachedSrcMapPtr_.resize(size_);
        cachedTgtAddress_.resize(size_);
        cachedTgtWeights_.resize(size_);
        cachedTgtWeightsSum_.resize(size_);
        cachedTgtMapPtr_.resize(size_);
    }
}


Foam::AMICache::AMICache(const bool toSource)
:
    size_(0),
    complete_(false),
    index0_(-1),
    index1_(-1),
    interpWeight_(0),
    coordSysPtr_(nullptr),
    theta_(),
    cachedSrcAddress_(),
    cachedSrcWeights_(),
    cachedSrcWeightsSum_(),
    cachedSrcMapPtr_(),
    cachedTgtAddress_(),
    cachedTgtWeights_(),
    cachedTgtWeightsSum_(),
    cachedTgtMapPtr_(),
    toSource_(toSource)
{}


Foam::AMICache::AMICache(const AMICache& cache)
:
    size_(cache.size_),
    complete_(cache.complete_),
    index0_(cache.index0_),
    index1_(cache.index1_),
    interpWeight_(cache.interpWeight_),
    coordSysPtr_(nullptr),  // Need to clone as cylindricalCS
    theta_(cache.theta_),
    cachedSrcAddress_(cache.cachedSrcAddress_),
    cachedSrcWeights_(cache.cachedSrcWeights_),
    cachedSrcWeightsSum_(cache.cachedSrcWeightsSum_),
    cachedSrcMapPtr_(cache.cachedSrcMapPtr_.size()),  // Need to clone
    cachedTgtAddress_(cache.cachedTgtAddress_),
    cachedTgtWeights_(cache.cachedTgtWeights_),
    cachedTgtWeightsSum_(cache.cachedTgtWeightsSum_),
    cachedTgtMapPtr_(cache.cachedTgtMapPtr_.size()),  // Need to clone
    toSource_(cache.toSource_)
{
    if (cache.coordSysPtr_)
    {
        coordSysPtr_.reset(new coordSystem::cylindrical(cache.coordSysPtr_()));
    }

    forAll(cachedSrcMapPtr_, cachei)
    {
        cachedSrcMapPtr_[cachei].reset(cache.cachedSrcMapPtr_[cachei].clone());
    }

    forAll(cachedTgtMapPtr_, cachei)
    {
        cachedTgtMapPtr_[cachei].reset(cache.cachedTgtMapPtr_[cachei].clone());
    }
}


Foam::AMICache::AMICache
(
    const AMICache& cache,
    const AMIInterpolation& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    size_(cache.size_),
    complete_(cache.complete_),
    index0_(cache.index0_),
    index1_(cache.index1_),
    interpWeight_(cache.interpWeight_),
    coordSysPtr_(nullptr),
    theta_(cache.theta_),
    cachedSrcAddress_(cache.size_),
    cachedSrcWeights_(cache.size_),
    cachedSrcWeightsSum_(cache.size_),
    cachedSrcMapPtr_(cache.size_),
    cachedTgtAddress_(cache.size_),
    cachedTgtWeights_(cache.size_),
    cachedTgtWeightsSum_(cache.size_),
    cachedTgtMapPtr_(cache.size_),
    toSource_(cache.toSource_)
{
    if (size_ > 0 && fineAMI.comm() != -1)
    {
        for (label cachei : {index0_, index1_})
        {
            if (cachei == -1) continue;

            scalarField dummySrcMagSf;
            labelListList srcAddress;
            scalarListList srcWeights;
            scalarField srcWeightsSum;
            autoPtr<mapDistribute> tgtMapPtr;

            AMIInterpolation::agglomerate
            (
                cache.cachedTgtMapPtr()[cachei],
                fineAMI.srcMagSf(),
                cache.cachedSrcAddress()[cachei],
                cache.cachedSrcWeights()[cachei],

                sourceRestrictAddressing,
                targetRestrictAddressing,

                dummySrcMagSf,
                srcAddress,
                srcWeights,
                srcWeightsSum,
                tgtMapPtr,
                fineAMI.comm()
            );

            scalarField dummyTgtMagSf;
            labelListList tgtAddress;
            scalarListList tgtWeights;
            scalarField tgtWeightsSum;
            autoPtr<mapDistribute> srcMapPtr;

            AMIInterpolation::agglomerate
            (
                cache.cachedSrcMapPtr()[cachei],
                fineAMI.tgtMagSf(),
                cache.cachedTgtAddress()[cachei],
                cache.cachedTgtWeights()[cachei],

                targetRestrictAddressing,
                sourceRestrictAddressing,

                dummyTgtMagSf,
                tgtAddress,
                tgtWeights,
                tgtWeightsSum,
                srcMapPtr,
                fineAMI.comm()
            );

            cachedSrcAddress_[cachei] = srcAddress;
            cachedSrcWeights_[cachei] = srcWeights;
            cachedSrcWeightsSum_[cachei] = srcWeightsSum;
            cachedSrcMapPtr_[cachei] = srcMapPtr.clone();

            cachedTgtAddress_[cachei] = tgtAddress;
            cachedTgtWeights_[cachei] = tgtWeights;
            cachedTgtWeightsSum_[cachei] = tgtWeightsSum;
            cachedTgtMapPtr_[cachei] = tgtMapPtr.clone();
        }
    }
}


Foam::AMICache::AMICache(Istream& is)
:
    size_(readLabel(is)),
    complete_(readBool(is)),

    index0_(-1),
    index1_(-1),
    interpWeight_(0),
    coordSysPtr_(nullptr),
    theta_(),
    cachedSrcAddress_(),
    cachedSrcWeights_(),
    cachedSrcWeightsSum_(),
    cachedSrcMapPtr_(),
    cachedTgtAddress_(),
    cachedTgtWeights_(),
    cachedTgtWeightsSum_(),
    cachedTgtMapPtr_()
{
    const bitSet goodMap(is);

    if (goodMap.size())
    {
        is >> index0_
           >> index1_
           >> interpWeight_
           >> theta_;

        const bool goodCoord(readBool(is));
        if (goodCoord)
        {
            coordSysPtr_.reset(new coordSystem::cylindrical(is));
        }

        is >> cachedSrcAddress_
           >> cachedSrcWeights_
           >> cachedSrcWeightsSum_;

        cachedSrcMapPtr_.setSize(goodMap.size());
        forAll(goodMap, cachei)
        {
            if (goodMap[cachei])
            {
                cachedSrcMapPtr_[cachei].reset(new mapDistribute(is));
            }
        }

        is >> cachedTgtAddress_
           >> cachedTgtWeights_
           >> cachedTgtWeightsSum_;

        cachedTgtMapPtr_.setSize(goodMap.size());
        forAll(goodMap, cachei)
        {
            if (goodMap[cachei])
            {
                cachedTgtMapPtr_[cachei].reset(new mapDistribute(is));
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::AMICache::addToCache
(
    const AMIInterpolation& ami,
    const point& globalPoint,
    const vector& rotationAxis,
    const vector& rotationCentre
)
{
    DebugPout<< "-- addToCache" << endl;


    if (!active())
    {
        DebugInfo<< "-- addToCache - deactivated" << endl;
        return;
    }

    if (!coordSysPtr_)
    {
        DebugInfo
            << "Creating rotation co-ordinate system:"
            << " rotationCentre:" << rotationCentre
            << " rotationAxis:" << rotationAxis
            << " p:" << globalPoint
            << endl;

        coordSysPtr_.reset
        (
            new coordSystem::cylindrical(rotationCentre, rotationAxis)
        );
        DebugPout<< "Coord sys:" << coordSysPtr_() << endl;
    }

    // Check if cache is complete
    if (!complete_)
    {
        for (const scalar angle : theta_)
        {
            if (angle > constant::mathematical::twoPi)
            {
                complete_ = false;
                break;
            }
        }
    }

    if (!complete_)
    {
        const scalar theta = getRotationAngle(globalPoint);

        const label bini = theta/constant::mathematical::twoPi*size_;

        DebugPout<< "--   bini:" << bini << " for theta:" << theta << endl;

        if (theta_[bini] > constant::mathematical::twoPi)
        {
            DebugPout<< "--   setting cache at index " << bini << endl;

            // New entry
            theta_[bini] = theta;

            cachedSrcAddress_[bini] = ami.srcAddress();
            cachedSrcWeights_[bini] = ami.srcWeights();
            cachedSrcWeightsSum_[bini] = ami.srcWeightsSum();

            if (ami.hasSrcMap())
            {
                cachedSrcMapPtr_[bini] = ami.srcMap().clone();
            }

            cachedTgtAddress_[bini] = ami.tgtAddress();
            cachedTgtWeights_[bini] = ami.tgtWeights();
            cachedTgtWeightsSum_[bini] = ami.tgtWeightsSum();

            if (ami.hasTgtMap())
            {
                cachedTgtMapPtr_[bini] = ami.tgtMap().clone();
            }
        }
    }
}


bool Foam::AMICache::restoreCache(const point& globalPoint)
{
    DebugPout<< "-- restoreCache" << endl;

    index0_ = -1;
    index1_ = -1;
    interpWeight_ = -1;

    if (!coordSysPtr_ || size_ == -1)
    {
        return false;
    }

    const scalar theta = getRotationAngle(globalPoint);
    const label bini = theta/constant::mathematical::twoPi*size_;

    DebugPout<< "--   bini:" << bini << " for theta:" << theta << endl;

    const auto validIndex = [&](const scalar bini)
    {
        return theta_[bini] < constant::mathematical::twoPi;
    };

    bool cacheValid = false;

    if (validIndex(bini))
    {
        // Find participating bins
        if (mag(theta - theta_[bini]) < cacheThetaTolerance_)
        {
            // Hit cached value - no interpolation needed
            index0_ = bini;
            cacheValid = true;
        }
        else if (theta > theta_[bini])
        {
            // Check that previous bin is valid
            const label i1 = theta_.fcIndex(bini);
            if (validIndex(i1))
            {
                index0_ = bini;
                index1_ = i1;
                cacheValid = true;
            }
        }
        else // (theta < theta_[bini])
        {
            // Check that previous bin is valid
            const label i1 = theta_.rcIndex(bini);
            if (validIndex(i1))
            {
                index0_ = i1;
                index1_ = bini;
                cacheValid = true;
            }
        }

        if (!cacheValid)
        {
            DebugPout<< "-- no cache available" << endl;
            return false;
        }


        // Calculate weighting factor
        if (index1_ != -1)
        {
            const scalar t0 = theta_[index0_];
            scalar t1 = theta_[index1_];

            if (index1_ < index0_)
            {
                t1 += constant::mathematical::twoPi;
            }

            // Set time-based weighting factor
            interpWeight_ = (theta - t0)/(t1 - t0);

            DebugInfo
                << "--   i0:" << index0_ << " i1:" << index1_
                << " w:" << interpWeight_ << endl;
        }
    }
    else
    {
        DebugPout<< "  -- no cache available" << endl;
    }

    return cacheValid;
}


void Foam::AMICache::write(Ostream& os) const
{
    if (size_ > 0)
    {
        os.writeEntry("cacheSize", size_);
    }
}


bool Foam::AMICache::writeData(Ostream& os) const
{
    os  << token::SPACE<< size_
        << token::SPACE<< complete_;

    bitSet goodMap(cachedSrcMapPtr_.size());
    forAll(goodMap, cachei)
    {
        goodMap.set(cachei, cachedSrcMapPtr_[cachei].good());
    }
    os  << token::SPACE << goodMap;

    if (goodMap.size())
    {
        os  << token::SPACE << index0_
            << token::SPACE << index1_
            << token::SPACE << interpWeight_
            << token::SPACE << theta_;

        os  << token::SPACE << coordSysPtr_.good();

        if (coordSysPtr_.good())
        {
            os  << token::SPACE << coordSysPtr_();
        }

        os  << token::SPACE << cachedSrcAddress_
            << token::SPACE << cachedSrcWeights_
            << token::SPACE << cachedSrcWeightsSum_;

        for (const auto& index : goodMap)
        {
            os  << token::SPACE << cachedSrcMapPtr_[index]();
        }

        os  << token::SPACE << cachedTgtAddress_
            << token::SPACE << cachedTgtWeights_
            << token::SPACE << cachedTgtWeightsSum_;

        for (const auto& index : goodMap)
        {
            os  << token::SPACE << cachedTgtMapPtr_[index]();
        }
    }

    return true;
}


// ************************************************************************* //
