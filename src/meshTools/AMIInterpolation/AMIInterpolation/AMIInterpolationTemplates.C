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
#include "AMIFieldOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolate");

    label inSize = cop.toSource() ? tgtAddress_.size() : srcAddress_.size();

    if (fld.size() != inSize)
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to expected field size ("
            << inSize << ")" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }


    label outSize = cop.toSource() ? srcAddress_.size() : tgtAddress_.size();
    result.resize_nocopy(outSize);

    if (distributed())
    {
        const auto& map = cop.toSource() ? tgtMapPtr_() : srcMapPtr_();
        List<Type> work = fld;  // deep copy
        map.distribute(work);

        // Apply interpolation
        cop(result, work, defaultValues);
    }
    else
    {
        // Apply interpolation
        cop(result, fld, defaultValues);
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolate
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New();

    interpolate(fld, cop, tresult.ref(), defaultValues);

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolate
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues,
    const lowWeightCorrectionBase::option& lwOption
) const
{
    AMICorrectedMultiplyWeightedOp<Type> cop(*this, true, lwOption);

    return interpolate(fld, cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues,
    const lowWeightCorrectionBase::option& lwOption
) const
{
    return interpolateToSource(tFld(), defaultValues, lwOption);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues,
    const lowWeightCorrectionBase::option& lwOption
) const
{
    AMICorrectedMultiplyWeightedOp<Type> cop(*this, false, lwOption);

    return interpolate(fld, cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues,
    const lowWeightCorrectionBase::option& lwOption
) const
{
    return interpolateToTarget(tFld(), defaultValues, lwOption);
}


// ************************************************************************* //
