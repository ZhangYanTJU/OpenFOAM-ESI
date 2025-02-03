/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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
    const bool toSource,
    const UList<Type>& fld,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    const auto& cAddress = (toSource ? cachedSrcAddress_ : cachedTgtAddress_);
    const auto& cWeights = (toSource ? cachedSrcWeights_ : cachedTgtWeights_);
    const auto& cWeightsSum =
    (
        toSource
      ? cachedSrcWeightsSum_
      : cachedTgtWeightsSum_
    );


    auto wsum = [&](List<Type>& res, const label i){
        weightedSum
        (
            lowWeightCorrection_,
            cAddress[i],
            cWeights[i],
            cWeightsSum[i],
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

        //result = (r1 - r0)*cachedWeight_ + r0;
        forAll(result, i)
        {
            result[i] = lerp(r0[i], r1[i], cachedWeight_);
        }
    }
    else
    {
        // Both -1 => equates to non-caching
        weightedSum
        (
            lowWeightCorrection_,
            (toSource ? srcAddress_ : tgtAddress_),
            (toSource ? srcWeights_ : tgtWeights_),
            (toSource ? srcWeightsSum_ : tgtWeightsSum_),
            fld,
            multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
            result,
            defaultValues
        );
    }
}


template<class Type, class CombineOp, class InterpolateOp>
void Foam::AMIInterpolation::interpolate
(
    const bool toSource,
    const UList<Type>& fld,
    const CombineOp& cop,
    const InterpolateOp& iop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    // Note: behaves as old AMIInterpolation::interpolateToSource if toSource=true

    // Get data locally and do a weighted sum

    addProfiling(ami, "AMIInterpolation::interpolate");

    auto checkSizes = [&](
        const UList<Type>& fld,
        const labelListList& srcAddr,
        const labelListList& tgtAddr,
        const UList<Type>& defVals
    )
    {
        const word srcName = toSource ? "source" : "target";
        const word tgtName = toSource ? "target" : "source";

        if (fld.size() != tgtAddr.size())
        {
            FatalErrorInFunction
                << "Supplied field size is not equal to "
                << tgtName << " patch size" << nl
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
                << " but number of default values is not equal to "
                << srcName << " patch size" << nl
                << "    default values = " << defVals.size() << nl
                << "    source patch   = " << srcAddr.size() << nl
                << abort(FatalError);
        }
    };


    // Work space for if distributed
    List<Type> work;


    List<Type> result0;
    if (cachedIndex0_ != -1)
    {
        result0 = result;

        const auto& srcAddress =
        (
            toSource
          ? cachedSrcAddress_[cachedIndex0_]
          : cachedTgtAddress_[cachedIndex0_]
        );
        const auto& srcWeights =
        (
            toSource
          ? cachedSrcWeights_[cachedIndex0_]
          : cachedTgtWeights_[cachedIndex0_]
        );
        const auto& srcWeightsSum =
        (
            toSource
          ? cachedSrcWeightsSum_[cachedIndex0_]
          : cachedTgtWeightsSum_[cachedIndex0_]
        );
        const auto& tgtAddress =
        (
            toSource
          ? cachedTgtAddress_[cachedIndex0_]
          : cachedSrcAddress_[cachedIndex0_]
        );

        checkSizes(fld, srcAddress, tgtAddress, defaultValues);

        if (distributed())
        {
            const mapDistribute& map =
            (
                toSource
              ? cachedTgtMapPtr_[cachedIndex0_]()
              : cachedSrcMapPtr_[cachedIndex0_]()
            );

            if (map.comm() == -1)
            {
                return;
            }
            
            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        result0.resize_nocopy(srcAddress.size());
        weightedSum
        (
            lowWeightCorrection_,
            srcAddress,
            srcWeights,
            srcWeightsSum,
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

        const auto& srcAddress =
        (
            toSource
          ? cachedSrcAddress_[cachedIndex1_]
          : cachedTgtAddress_[cachedIndex1_]
        );
        const auto& srcWeights =
        (
            toSource
          ? cachedSrcWeights_[cachedIndex1_]
          : cachedTgtWeights_[cachedIndex1_]
        );
        const auto& srcWeightsSum =
        (
            toSource
          ? cachedSrcWeightsSum_[cachedIndex1_]
          : cachedTgtWeightsSum_[cachedIndex1_]
        );
        const auto& tgtAddress =
        (
            toSource
          ? cachedTgtAddress_[cachedIndex1_]
          : cachedSrcAddress_[cachedIndex1_]
        );

        checkSizes(fld, srcAddress, tgtAddress, defaultValues);

        if (distributed())
        {
            const mapDistribute& map =
            (
                toSource
              ? cachedTgtMapPtr_[cachedIndex1_]()
              : cachedSrcMapPtr_[cachedIndex1_]()
            );

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        result1.resize_nocopy(srcAddress.size());
        weightedSum
        (
            lowWeightCorrection_,
            srcAddress,
            srcWeights,
            srcWeightsSum,
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
        forAll(result, i)
        {
            iop(result[i], i, i, result0[i], i, result1[i], cachedWeight_);
        }
    }
    else
    {
        // No cache - evaluate the AMI

        const auto& srcAddress = (toSource ? srcAddress_ : tgtAddress_);
        const auto& srcWeights = (toSource ? srcWeights_ : tgtWeights_);
        const auto& srcWeightsSum =
            (toSource ? srcWeightsSum_ : tgtWeightsSum_);
        const auto& tgtAddress = (toSource ? tgtAddress_ : srcAddress_);

        checkSizes(fld, srcAddress, tgtAddress, defaultValues);

        if (distributed() && tgtMapPtr_)
        {
            const mapDistribute& map =
            (
                toSource
              ? tgtMapPtr_()
              : srcMapPtr_()
            );

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        result.resize_nocopy(srcAddress.size());
        weightedSum
        (
            lowWeightCorrection_,
            srcAddress,
            srcWeights,
            srcWeightsSum,
            (distributed() ? work : fld),
            cop,
            result,
            defaultValues
        );
    }
}


// Leave API intact below!
template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    // In-place interpolation

    addProfiling(ami, "AMIInterpolation::interpolateToTarget");

    // Wrap lerp operator to operate inplace
    auto iop = [&]
    (
        Type& res,
        const label i,
        const label ia,
        const Type& a,
        const label ib,
        const Type& b,
        const scalar w
    )
    {
        res = lerp(a, b, w);
    };

    interpolate
    (
        false,                  // interpolate to target
        fld,
        cop,
        iop,
        result,
        defaultValues
    );
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
    // In-place interpolation

    addProfiling(ami, "AMIInterpolation::interpolateToSource");

    // Wrap lerp operator to operate inplace
    auto iop = [&]
    (
        Type& res,
        const label i,
        const label ia,
        const Type& a,
        const label ib,
        const Type& b,
        const scalar w
    )
    {
        res = lerp(a, b, w);
    };


    interpolate
    (
        true,                   // toSource,
        fld,
        cop,
        iop,
        result,
        defaultValues
    );
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
