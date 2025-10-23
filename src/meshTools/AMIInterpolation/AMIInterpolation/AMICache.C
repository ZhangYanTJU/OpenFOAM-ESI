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
    rotationAxis_(dict.getOrDefault<vector>("rotationAxis", Zero)),
    rotationCentre_(dict.getOrDefault<point>("rotationCentre", Zero)),
    maxThetaDeg_(dict.getOrDefault<scalar>("maxThetaDeg", 5)),
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
    rotationAxis_(Zero),
    rotationCentre_(Zero),
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
    rotationAxis_(cache.rotationAxis_),
    rotationCentre_(cache.rotationCentre_),
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
    rotationAxis_(cache.rotationAxis_),
    rotationCentre_(cache.rotationCentre_),
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

    rotationAxis_(is),
    rotationCentre_(is),
    maxThetaDeg_(readScalar(is)),

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
    const point& globalPoint
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
            << " rotationCentre:" << rotationCentre_
            << " rotationAxis:" << rotationAxis_
            << " p:" << globalPoint
            << endl;

        coordSysPtr_.reset
        (
            new coordSystem::cylindrical(rotationCentre_, rotationAxis_)
        );
        DebugPout<< "Coord sys:" << coordSysPtr_() << endl;
    }

    // Check if cache is complete
    if (!complete_)
    {
        for (const scalar angle : theta_)
        {
            // Check against null value; if any are present, cache is incomplete
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

        DebugPout<< "  -- bini:" << bini << " for theta:" << theta << endl;

        // Check if already have entry for this bin
        if (theta_[bini] > constant::mathematical::twoPi)
        {
            DebugPout<< "  -- setting cache at index " << bini << endl;

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
    const scalar twoPi = constant::mathematical::twoPi;
    const label bini = theta/twoPi*size_;

    DebugPout<< "  -- bini:" << bini << " for theta:" << theta << endl;

    const auto validIndex = [&](const scalar bini)
    {
        return theta_[bini] < constant::mathematical::twoPi;
    };

    // Maximum angle in degrees for which to search for cached bins
    const scalar maxThetaStencil = maxThetaDeg_*constant::mathematical::pi/180.0;

    if
    (
        validIndex(bini)
     && (
            mag(theta - theta_[bini]) < cacheThetaTolerance_
         || mag(theta - twoPi - theta_[bini]) < cacheThetaTolerance_
        )
    )
    {
        // Hit cached value - no interpolation needed
        // index1_ = -1 indicates no interpolation
        index0_ = bini;
        index1_ = -1;
        interpWeight_ = 0;

        DebugInfo
            << "  -- t0:" << theta_[index0_] << " theta:" << theta
            << " i0:" << index0_ << " i1:" << index1_
            << " w:" << interpWeight_ << endl;
        return true;
    }
    else
    {
        // Find indices and values bracketing theta
        const label nBin = theta_.size();

        // Participating theta values and bin addresses
        // - Note we add wrap-around values at start and end
        DynamicList<scalar> thetap(nBin+2);
        DynamicList<label> binAddresses(nBin+2);

        // Initialise wrap-around values
        thetap.push_back(0);
        binAddresses.push_back(-1);
        forAll(theta_, thetai)
        {
            if (validIndex(thetai))
            {
                thetap.push_back(theta_[thetai]);
                binAddresses.push_back(thetai);
            }
        }

        // Check that we have enough data points for interpolation
        // - We added storage for lower wrap-around value, and we then need
        //   at least 2 additional values for the interpolation
        if (thetap.size() < 3)
        {
            DebugPout<< "  -- no cache available" << endl;
            return false;
        }

        // Set wrap-around values if we have sufficient data
        thetap[0] = thetap.last() - twoPi;
        binAddresses[0] = binAddresses.last();
        thetap.push_back(thetap[1] + twoPi);
        binAddresses.push_back(binAddresses[1]);

        // Find interpolation indices
        label loweri = labelMax;
        label i = 0;
        while (i < thetap.size())
        {
            if (thetap[i] < theta)
            {
                loweri = i;
            }
            else
            {
                break;
            }
            ++i;
        }

        if (loweri == labelMax)
        {
            DebugPout<< "  -- no cache available" << endl;
            return false;
        }

        label upperi = labelMax;
        i = thetap.size() - 1;
        while (i >= 0)
        {
            if (thetap[i] > theta)
            {
                upperi = i;
            }
            else
            {
                break;
            }
            --i;
        }

        if (upperi == labelMax)
        {
            DebugPout<< "  -- no cache available" << endl;
            return false;
        }

        // Ensure distances are valid
        if (upperi == loweri)
        {
            DebugPout
                << "  -- no cache available: theta:" << theta
                << " lower:" << loweri << " upper:" << upperi << endl;
            return false;
        }
        if (mag(theta - thetap[loweri]) > maxThetaStencil)
        {
            DebugPout
                << "  -- no cache available: theta:" << theta
                << " lower:" << thetap[loweri] << endl;
            return false;
        }
        if (mag(theta - thetap[upperi]) > maxThetaStencil)
        {
            DebugPout
                << "  -- no cache available: theta:" << theta
                << " upper:" << thetap[upperi] << endl;
            return false;
        }

        index0_ = binAddresses[loweri];
        index1_ = binAddresses[upperi];
        interpWeight_ =
            (theta - theta_[index0_])/(theta_[index1_] - theta_[index0_]);

        DebugInfo
            << theta_.size()
            << "  -- t0:" << theta_[index0_] << " theta:" << theta
            << " t1:" << theta_[index1_]
            << " i0:" << index0_ << " i1:" << index1_
            << " w:" << interpWeight_ << endl;


        return true;
    }

    // If we get here then no valid cache found within stencil
    DebugPout<< "  -- no cache available" << endl;
    return false;
}


void Foam::AMICache::write(Ostream& os) const
{
    if (size_ > 0)
    {
        os.writeEntry("cacheSize", size_);
        os.writeEntry("rotationAxis", rotationAxis_);
        os.writeEntry("rotationCentre", rotationCentre_);
        os.writeEntry("maxThetaDeg", maxThetaDeg_);
    }
}


bool Foam::AMICache::writeData(Ostream& os) const
{
    os  << token::SPACE<< size_
        << token::SPACE<< rotationAxis_
        << token::SPACE<< rotationCentre_
        << token::SPACE<< maxThetaDeg_
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
