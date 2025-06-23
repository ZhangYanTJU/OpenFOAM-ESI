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
    // Note: using non-caching AMI
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

    cache_.setDirection(toSource);

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
    if (cache_.index0() != -1)
    {
        result0 = result;

        const auto& srcAddress = cache_.cSrcAddress0();
        const auto& srcWeights = cache_.cSrcWeights0();
        const auto& srcWeightsSum = cache_.cSrcWeightsSum0();
        const auto& tgtAddress = cache_.cTgtAddress0();

        checkSizes(fld, srcAddress, tgtAddress, defaultValues);

        if (distributed() && cache_.cTgtMapPtr0())
        {
            const mapDistribute& map = cache_.cTgtMapPtr0()();

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        if constexpr (is_contiguous_scalar<Type>::value)
        {
            result0 = Zero;
        }

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
    if (cache_.index1() != -1)
    {
        result1 = result;

        const auto& srcAddress = cache_.cSrcAddress1();
        const auto& srcWeights = cache_.cSrcWeights1();
        const auto& srcWeightsSum = cache_.cSrcWeightsSum1();
        const auto& tgtAddress = cache_.cTgtAddress1();

        checkSizes(fld, srcAddress, tgtAddress, defaultValues);

        if (distributed() && cache_.cTgtMapPtr1())
        {
            const mapDistribute& map = cache_.cTgtMapPtr1()();

            if (map.comm() == -1)
            {
                return;
            }

            work.resize_nocopy(map.constructSize());
            SubList<Type>(work, fld.size()) = fld;  // deep copy
            map.distribute(work);
        }

        if constexpr (is_contiguous_scalar<Type>::value)
        {
            result1 = Zero;
        }

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

    if (cache_.applyLower())
    {
        result = result0;
    }
    else if (cache_.applyUpper())
    {
        result = result1;
    }
    else if (cache_.applyInterpolate())
    {
        forAll(result, i)
        {
            iop(result[i], i, i, result0[i], i, result1[i], cache_.weight());
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

        if constexpr (is_contiguous_scalar<Type>::value)
        {
            result = Zero;
        }

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
