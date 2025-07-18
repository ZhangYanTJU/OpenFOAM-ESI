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

    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }
    else if
    (
        (lowWeightCorrection_ > 0)
     && (defaultValues.size() != tgtAddress_.size())
    )
    {
        FatalErrorInFunction
            << "Employing default values when sum of weights falls below "
            << lowWeightCorrection_
            << " but supplied default field size is not equal to target "
            << "patch size" << nl
            << "    default values = " << defaultValues.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << abort(FatalError);
    }

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

    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }
    else if
    (
        (lowWeightCorrection_ > 0)
     && (defaultValues.size() != srcAddress_.size())
    )
    {
        FatalErrorInFunction
            << "Employing default values when sum of weights falls below "
            << lowWeightCorrection_
            << " but number of default values is not equal to source "
            << "patch size" << nl
            << "    default values = " << defaultValues.size() << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << abort(FatalError);
    }

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
