/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "profiling.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolation::weightedSum
(
    const scalar lowWeightCorrection,
    const labelListList& allSlots,
    const scalarListList& allWeights,
    const scalarField& weightsSum,
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
)
{
//     DebugVar("AMIInterpolation::weightedSum");

// Info<< "allSlots.size():" << allSlots.size() << nl
//     << "allWeights.size():" << allWeights.size() << nl
//     << "weightsSum.size():" << weightsSum.size() << nl
//     << "fld.size():" << fld.size() << nl
//     << "defaultValues.size():" << defaultValues.size() << nl;

    if (lowWeightCorrection > 0)
    {
        forAll(result, facei)
        {
            if (weightsSum[facei] < lowWeightCorrection)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& slots = allSlots[facei];
                const scalarList& weights = allWeights[facei];

                forAll(slots, i)
                {
                    cop(result[facei], facei, fld[slots[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& slots = allSlots[facei];
            const scalarList& weights = allWeights[facei];

            forAll(slots, i)
            {
                cop(result[facei], facei, fld[slots[i]], weights[i]);
            }
        }
    }
}


template<class Type>
void Foam::AMIInterpolation::weightedSum
(
    const bool interpolateToSource,
    const UList<Type>& fld,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    // DebugVar("AMIInterpolation::weightedSum");
// Info<< "cachedIndex0:" << cachedIndex0_ << " cachedIndex1:" << cachedIndex1_ << endl;

// Info<< "cachedSrcAddress_.size():" << cachedSrcAddress_.size() << nl
//     << "cachedTgtAddress_.size():" << cachedTgtAddress_.size() << nl
//     << "cachedSrcWeights_.size():" << cachedSrcWeights_.size() << nl
//     << "cachedTgtWeights_.size():" << cachedTgtWeights_.size() << nl
//     << "cachedSrcWeightsSum_.size():" << cachedSrcWeightsSum_.size() << nl
//     << "cachedTgtWeightsSum_.size():" << cachedTgtWeightsSum_.size() << nl;

    auto wsum = [&](List<Type>& res, const label i){
        weightedSum
        (
            lowWeightCorrection_,
            (interpolateToSource ? cachedSrcAddress_[i] : cachedTgtAddress_[i]),
            (interpolateToSource ? cachedSrcWeights_[i] : cachedTgtWeights_[i]),
            (interpolateToSource ? scalarField(cachedSrcWeightsSum_[i]) : scalarField(cachedTgtWeightsSum_[i])),
            fld,
            multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
            res,
            defaultValues
        );
    };

    if (cachedIndex0_ != -1 && cachedIndex1_ == -1)
    {
        wsum(result, cachedIndex0_);
    }
    else if (cachedIndex0_ == -1 && cachedIndex1_ != -1)
    {
        wsum(result, cachedIndex1_);
    }
    else if (cachedIndex0_ != -1 && cachedIndex1_ != -1)
    {
        List<Type> r0(result);
        wsum(r0, cachedIndex0_);
        List<Type> r1(result);
        wsum(r1, cachedIndex1_);
        result = (r1 - r0)*cachedWeight_ + r0;
    }
    else
    {
        // Both -1 => equates to non-caching
        weightedSum
        (
            lowWeightCorrection_,
            (interpolateToSource ? srcAddress_ : tgtAddress_),
            (interpolateToSource ? srcWeights_ : tgtWeights_),
            (interpolateToSource ? srcWeightsSum_ : tgtWeightsSum_),
            fld,
            multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToTarget");

    auto checkSizes = [&](
        const UList<Type>& fld,
        const labelListList& srcAddr,
        const labelListList& tgtAddr,
        const UList<Type>& defVals
    )
    {
        if (fld.size() != srcAddr.size())
        {
            FatalErrorInFunction
                << "Supplied field size is not equal to source patch size" << nl
                << "    source patch   = " << srcAddr.size() << nl
                << "    target patch   = " << tgtAddr.size() << nl
                << "    supplied field = " << fld.size()
                << abort(FatalError);
        }
        else if
        (
            (lowWeightCorrection_ > 0) && (defVals.size() != tgtAddr.size())
        )
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but number of default values is not equal to target "
                << "patch size" << nl
                << "    default values = " << defVals.size() << nl
                << "    target patch   = " << tgtAddr.size() << nl
                << abort(FatalError);
        }
    };

    List<Type> result0;
    if (cachedIndex0_ != -1)
    {
        result0 = result;

        const auto& srcAddr = cachedSrcAddress_[cachedIndex0_];
        const auto& tgtAddr = cachedTgtAddress_[cachedIndex0_];

        checkSizes(fld, srcAddr, tgtAddr, defaultValues);

        result0.setSize(tgtAddr.size());
        List<Type> work;

        if (distributed() && cachedSrcMapPtr_[cachedIndex0_])
        {
            const mapDistribute& map = cachedSrcMapPtr_[cachedIndex0_];

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            tgtAddr,
            cachedTgtWeights_[cachedIndex0_],
            scalarField(cachedTgtWeightsSum_[cachedIndex0_]),
            (distributed() ? work : fld),
            cop,
            result0,
            defaultValues
        );
    }

    List<Type> result1;
    if (cachedIndex1_ != -1)
    {
        result1 = result;

        const auto& srcAddr = cachedSrcAddress_[cachedIndex1_];
        const auto& tgtAddr = cachedTgtAddress_[cachedIndex1_];

        checkSizes(fld, srcAddr, tgtAddr, defaultValues);

        result1.setSize(tgtAddr.size());
        List<Type> work;

        if (distributed() && cachedSrcMapPtr_[cachedIndex1_])
        {
            const mapDistribute& map = cachedSrcMapPtr_[cachedIndex1_];

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            tgtAddr,
            cachedTgtWeights_[cachedIndex1_],
            scalarField(cachedTgtWeightsSum_[cachedIndex1_]),
            (distributed() ? work : fld),
            cop,
            result1,
            defaultValues
        );
    }

    if (cachedIndex0_ != -1 && cachedIndex1_ == -1)
    {
        result = result0;
    }
    else if (cachedIndex0_ == -1 && cachedIndex1_ != -1)
    {
        result = result1;
    }
    else if (cachedIndex0_ != -1 && cachedIndex1_ != -1)
    {
        result = (result1 - result0)*cachedWeight_ + result0;
    }
    else
    {
        // No cache - evaluate the AMI
        checkSizes(fld, srcAddress_, tgtAddress_, defaultValues);

        result.setSize(tgtAddress_.size());
        List<Type> work;

        if (distributed() && srcMapPtr_)
        {
            const mapDistribute& map = srcMapPtr_();

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            tgtAddress_,
            tgtWeights_,
            tgtWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToSource");

    auto checkSizes = [&](
        const UList<Type>& fld,
        const labelListList& srcAddr,
        const labelListList& tgtAddr,
        const UList<Type>& defVals
    )
    {
        if (fld.size() != tgtAddr.size())
        {
            FatalErrorInFunction
                << "Supplied field size is not equal to target patch size" << nl
                << "    source patch   = " << srcAddr.size() << nl
                << "    target patch   = " << tgtAddr.size() << nl
                << "    supplied field = " << fld.size()
                << abort(FatalError);
        }
        else if
        (
            (lowWeightCorrection_ > 0) && (defVals.size() != srcAddr.size())
        )
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but number of default values is not equal to source "
                << "patch size" << nl
                << "    default values = " << defVals.size() << nl
                << "    source patch   = " << srcAddr.size() << nl
                << abort(FatalError);
        }
    };

    List<Type> result0;
    if (cachedIndex0_ != -1)
    {
        result0 = result;

        const auto& srcAddr = cachedSrcAddress_[cachedIndex0_];
        const auto& tgtAddr = cachedTgtAddress_[cachedIndex0_];

        checkSizes(fld, srcAddr, tgtAddr, defaultValues);

        result0.setSize(srcAddr.size());
        List<Type> work;

        if (distributed() && cachedTgtMapPtr_[cachedIndex0_])
        {
            const mapDistribute& map = cachedTgtMapPtr_[cachedIndex0_];

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            srcAddr,
            cachedSrcWeights_[cachedIndex0_],
            scalarField(cachedSrcWeightsSum_[cachedIndex0_]),
            (distributed() ? work : fld),
            cop,
            result0,
            defaultValues
        );
    }

    List<Type> result1;
    if (cachedIndex1_ != -1)
    {
        result1 = result;

        const auto& srcAddr = cachedSrcAddress_[cachedIndex1_];
        const auto& tgtAddr = cachedTgtAddress_[cachedIndex1_];

        checkSizes(fld, srcAddr, tgtAddr, defaultValues);

        result0.setSize(srcAddr.size());
        List<Type> work;

        if (distributed() && cachedTgtMapPtr_[cachedIndex1_])
        {
            const mapDistribute& map = cachedTgtMapPtr_[cachedIndex1_];

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            srcAddr,
            cachedSrcWeights_[cachedIndex1_],
            scalarField(cachedSrcWeightsSum_[cachedIndex1_]),
            (distributed() ? work : fld),
            cop,
            result1,
            defaultValues
        );
    }

    if (cachedIndex0_ != -1 && cachedIndex1_ == -1)
    {
        result = result0;
    }
    else if (cachedIndex0_ == -1 && cachedIndex1_ != -1)
    {
        result = result1;
    }
    else if (cachedIndex0_ != -1 && cachedIndex1_ != -1)
    {
        result = (result1 - result0)*cachedWeight_ + result0;
    }
    else
    {
        // No cache - evaluate the AMI
        checkSizes(fld, srcAddress_, tgtAddress_, defaultValues);

        result.setSize(srcAddress_.size());
        List<Type> work;

        if (distributed() && tgtMapPtr_)
        {
            const mapDistribute& map = tgtMapPtr_();

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        weightedSum
        (
            lowWeightCorrection_,
            srcAddress_,
            srcWeights_,
            srcWeightsSum_,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(srcAddress_.size(), Zero);

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(tgtAddress_.size(), Zero);

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


// ************************************************************************* //
