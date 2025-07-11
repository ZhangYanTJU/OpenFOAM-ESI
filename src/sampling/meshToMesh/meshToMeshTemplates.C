/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "volFields.H"
#include "directFvPatchFieldMapper.H"
#include "calculatedFvPatchField.H"
#include "fvcGrad.H"
#include "distributedWeightedFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::meshToMesh::add
(
    UList<Type>& fld,
    const label offset
) const
{
    forAll(fld, i)
    {
        fld[i] += offset;
    }
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapSrcToTgt
(
    const UList<Type>& srcField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != tgtToSrcCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (distributed())
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        forAll(result, celli)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[celli];
            const scalarList& srcWeight = tgtToSrcCellWght_[celli];

            if (srcAddress.size())
            {
//                result[celli] = Zero;
                result[celli] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[celli], celli, work[srcI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, celli)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[celli];
            const scalarList& srcWeight = tgtToSrcCellWght_[celli];

            if (srcAddress.size())
            {
//                result[celli] = Zero;
                result[celli] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[celli], celli, srcField[srcI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapSrcToTgt
(
    const UList<Type>& srcField,
    const UList<typename outerProduct<vector, Type>::type>& srcGradField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != tgtToSrcCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (distributed())
    {
        if (returnReduceAnd(tgtToSrcCellVec_.empty()))
        {
            // No correction vectors calculated. Fall back to first order.
            mapSrcToTgt(srcField, cop, result);
            return;
        }

        const mapDistribute& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        List<typename outerProduct<vector, Type>::type> workGrad
        (
            srcGradField
        );
        map.distribute(workGrad);

        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];
            const pointList& srcVec = tgtToSrcCellVec_[cellI];

            if (srcAddress.size())
            {
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    const vector& v = srcVec[i];
                    const Type srcVal = work[srcI]+(workGrad[srcI]&v);
                    cbop(result[cellI], cellI, srcVal, w);
                }
            }
        }
    }
    else
    {
        if (tgtToSrcCellVec_.empty())
        {
            // No correction vectors calculated. Fall back to first order.
            mapSrcToTgt(srcField, cop, result);
            return;
        }

        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];
            const pointList& srcVec = tgtToSrcCellVec_[cellI];

            if (srcAddress.size())
            {
                // Do non-conservative interpolation
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    const vector& v = srcVec[i];
                    const Type srcVal = srcField[srcI]+(srcGradField[srcI]&v);
                    cbop(result[cellI], cellI, srcVal, w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField,
    const CombineOp& cop
) const
{
    auto tresult = tmp<Field<Type>>::New(tgtToSrcCellAddr_.size(), Zero);

    mapSrcToTgt(srcField, cop, tresult.ref());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const tmp<Field<Type>>& tsrcField,
    const CombineOp& cop
) const
{
    return mapSrcToTgt(tsrcField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField
) const
{
    return mapSrcToTgt(srcField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const tmp<Field<Type>>& tsrcField
) const
{
    return mapSrcToTgt(tsrcField());
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapTgtToSrc
(
    const UList<Type>& tgtField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != srcToTgtCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (distributed())
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        forAll(result, celli)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[celli];
            const scalarList& tgtWeight = srcToTgtCellWght_[celli];

            if (tgtAddress.size())
            {
                result[celli] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[celli], celli, work[tgtI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, celli)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[celli];
            const scalarList& tgtWeight = srcToTgtCellWght_[celli];

            if (tgtAddress.size())
            {
                result[celli] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[celli], celli, tgtField[tgtI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapTgtToSrc
(
    const UList<Type>& tgtField,
    const UList<typename outerProduct<vector, Type>::type>& tgtGradField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != srcToTgtCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (distributed())
    {
        if (returnReduceAnd(srcToTgtCellVec_.empty()))
        {
            // No correction vectors calculated. Fall back to first order.
            mapTgtToSrc(tgtField, cop, result);
            return;
        }

        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        List<typename outerProduct<vector, Type>::type> workGrad
        (
            tgtGradField
        );
        map.distribute(workGrad);

        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];
            const pointList& tgtVec = srcToTgtCellVec_[cellI];

            if (tgtAddress.size())
            {
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    const vector& v = tgtVec[i];
                    const Type tgtVal = work[tgtI]+(workGrad[tgtI]&v);
                    cbop(result[cellI], cellI, tgtVal, w);
                }
            }
        }
    }
    else
    {
        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];
            const pointList& tgtVec = srcToTgtCellVec_[cellI];

            if (tgtAddress.size())
            {
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    const vector& v = tgtVec[i];
                    const Type tgtVal = tgtField[tgtI]+(tgtGradField[tgtI]&v);
                    cbop(result[cellI], cellI, tgtVal, w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField,
    const CombineOp& cop
) const
{
    auto tresult = tmp<Field<Type>>::New(srcToTgtCellAddr_.size(), Zero);

    mapTgtToSrc(tgtField, cop, tresult.ref());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const tmp<Field<Type>>& ttgtField,
    const CombineOp& cop
) const
{
    return mapTgtToSrc(ttgtField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField
) const
{
    return mapTgtToSrc(tgtField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const tmp<Field<Type>>& ttgtField
) const
{
    return mapTgtToSrc(ttgtField(), plusEqOp<Type>());
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapInternalSrcToTgt
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    VolumeField<Type>& result,
    const bool secondOrder
) const
{
    if (secondOrder && returnReduceOr(tgtToSrcCellVec_.size()))
    {
        mapSrcToTgt
        (
            field,
            fvc::grad(field)().primitiveField(),
            cop,
            result.primitiveFieldRef()
        );
    }
    else
    {
        mapSrcToTgt(field, cop, result.primitiveFieldRef());
    }
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapAndOpSrcToTgt
(
    const AMIPatchToPatchInterpolation& AMI,
    const Field<Type>& srcField,
    Field<Type>& tgtField,
    const CombineOp& cop
) const
{
    tgtField = Type(Zero);

    AMI.interpolateToTarget
    (
        srcField,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tgtField,
        UList<Type>::null()
    );
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapSrcToTgt
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    VolumeField<Type>& result,
    const bool secondOrder
) const
{
    mapInternalSrcToTgt(field, cop, result, secondOrder);

    const PtrList<AMIPatchToPatchInterpolation>& AMIList = patchAMIs();

    auto& resultBf = result.boundaryFieldRef();

    forAll(AMIList, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        const fvPatchField<Type>& srcField = field.boundaryField()[srcPatchi];
        fvPatchField<Type>& tgtField = resultBf[tgtPatchi];

        // Clone and map (since rmap does not do general mapping)
        tmp<fvPatchField<Type>> tnewTgt
        (
            fvPatchField<Type>::New
            (
                srcField,
                tgtField.patch(),
                result(),
                distributedWeightedFvPatchFieldMapper
                (
                    AMIList[i].singlePatchProc(),
                    AMIList[i].comm(),      // for testing only
                    AMIList[i].hasSrcMap(), // pointer to map
                    AMIList[i].tgtAddress(),
                    AMIList[i].tgtWeights()
                )
            )
        );

        // Transfer all mapped quantities (value and e.g. gradient) onto
        // tgtField. Value will get overwritten below.
        tgtField.rmap(tnewTgt(), identity(tgtField.size()));

        // Override value to account for CombineOp (note: is dummy template
        // specialisation for plusEqOp)
        mapAndOpSrcToTgt(AMIList[i], srcField, tgtField, cop);
    }

    forAll(cuttingPatches_, i)
    {
        label patchi = cuttingPatches_[i];
        fvPatchField<Type>& pf = resultBf[patchi];
        pf == pf.patchInternalField();
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapSrcToTgt
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    const bool secondOrder
) const
{
    const fvMesh& tgtMesh = static_cast<const fvMesh&>(tgtRegion_);

    const fvBoundaryMesh& tgtBm = tgtMesh.boundary();
    const auto& srcBfld = field.boundaryField();

    PtrList<fvPatchField<Type>> tgtPatchFields(tgtBm.size());

    // construct tgt boundary patch types as copy of 'field' boundary types
    // note: this will provide place holders for fields with additional
    // entries, but these values will need to be reset
    forAll(tgtPatchID_, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        if (!tgtPatchFields.set(tgtPatchi))
        {
            tgtPatchFields.set
            (
                tgtPatchi,
                fvPatchField<Type>::New
                (
                    srcBfld[srcPatchi],
                    tgtMesh.boundary()[tgtPatchi],
                    fvPatchField<Type>::Internal::null(),
                    directFvPatchFieldMapper
                    (
                        labelList(tgtMesh.boundary()[tgtPatchi].size(), -1)
                    )
                )
            );
        }
    }

    // Any unset tgtPatchFields become calculated
    forAll(tgtPatchFields, tgtPatchi)
    {
        if (!tgtPatchFields.set(tgtPatchi))
        {
            // Note: use factory New method instead of direct generation of
            //       calculated so we keep constraints
            tgtPatchFields.set
            (
                tgtPatchi,
                fvPatchField<Type>::New
                (
                    fvPatchFieldBase::calculatedType(),
                    tgtMesh.boundary()[tgtPatchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
    }

    auto tresult =
        tmp<VolumeField<Type>>::New
        (
            tgtMesh.newIOobject
            (
                IOobject::scopedName
                (
                    type(),
                    "interpolate(" + field.name() + ")"
                )
            ),
            tgtMesh,
            field.dimensions(),
            Field<Type>(tgtMesh.nCells(), Zero),
            tgtPatchFields
        );

    mapSrcToTgt(field, cop, tresult.ref(), secondOrder);

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapSrcToTgt
(
    const tmp<VolumeField<Type>>& tfield,
    const CombineOp& cop,
    const bool secondOrder
) const
{
    return mapSrcToTgt(tfield(), cop, secondOrder);
}


template<class Type>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapSrcToTgt
(
    const VolumeField<Type>& field,
    const bool secondOrder
) const
{
    return mapSrcToTgt(field, plusEqOp<Type>(), secondOrder);
}


template<class Type>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapSrcToTgt
(
    const tmp<VolumeField<Type>>& tfield,
    const bool secondOrder
) const
{
    return mapSrcToTgt(tfield(), plusEqOp<Type>(), secondOrder);
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapInternalTgtToSrc
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    VolumeField<Type>& result,
    const bool secondOrder
) const
{
    if (secondOrder && returnReduceOr(srcToTgtCellVec_.size()))
    {
        mapTgtToSrc
        (
            field,
            fvc::grad(field)().primitiveField(),
            cop,
            result.primitiveFieldRef()
        );
    }
    else
    {
        mapTgtToSrc(field, cop, result.primitiveFieldRef());
    }
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapAndOpTgtToSrc
(
    const AMIPatchToPatchInterpolation& AMI,
    Field<Type>& srcField,
    const Field<Type>& tgtField,
    const CombineOp& cop
) const
{
    srcField = Type(Zero);

    AMI.interpolateToSource
    (
        tgtField,
        multiplyWeightedOp<Type, CombineOp>(cop),
        srcField,
        UList<Type>::null()
    );
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapTgtToSrc
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    VolumeField<Type>& result,
    const bool secondOrder
) const
{
    mapInternalTgtToSrc(field, cop, result, secondOrder);

    const PtrList<AMIPatchToPatchInterpolation>& AMIList = patchAMIs();

    forAll(AMIList, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        fvPatchField<Type>& srcField = result.boundaryFieldRef()[srcPatchi];
        const fvPatchField<Type>& tgtField = field.boundaryField()[tgtPatchi];

        // Clone and map (since rmap does not do general mapping)
        tmp<fvPatchField<Type>> tnewSrc
        (
            fvPatchField<Type>::New
            (
                tgtField,
                srcField.patch(),
                result(),
                distributedWeightedFvPatchFieldMapper
                (
                    AMIList[i].singlePatchProc(),
                    AMIList[i].comm(),      // only used for testing
                    AMIList[i].hasTgtMap(), // pointer to map
                    AMIList[i].srcAddress(),
                    AMIList[i].srcWeights()
                )
            )
        );
        // Transfer all mapped quantities (value and e.g. gradient) onto
        // srcField. Value will get overwritten below
        srcField.rmap(tnewSrc(), identity(srcField.size()));

        // Override value to account for CombineOp (could be dummy for
        // plusEqOp)
        mapAndOpTgtToSrc(AMIList[i], srcField, tgtField, cop);
    }

    forAll(cuttingPatches_, i)
    {
        label patchi = cuttingPatches_[i];
        fvPatchField<Type>& pf = result.boundaryFieldRef()[patchi];
        pf == pf.patchInternalField();
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapTgtToSrc
(
    const VolumeField<Type>& field,
    const CombineOp& cop,
    const bool secondOrder
) const
{
    const fvMesh& srcMesh = static_cast<const fvMesh&>(srcRegion_);

    const fvBoundaryMesh& srcBm = srcMesh.boundary();
    const auto& tgtBfld = field.boundaryField();

    PtrList<fvPatchField<Type>> srcPatchFields(srcBm.size());

    // construct src boundary patch types as copy of 'field' boundary types
    // note: this will provide place holders for fields with additional
    // entries, but these values will need to be reset
    forAll(srcPatchID_, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        if (!srcPatchFields.set(srcPatchi))
        {
            srcPatchFields.set
            (
                srcPatchi,
                fvPatchField<Type>::New
                (
                    tgtBfld[tgtPatchi],
                    srcMesh.boundary()[srcPatchi],
                    fvPatchField<Type>::Internal::null(),
                    directFvPatchFieldMapper
                    (
                        labelList(srcMesh.boundary()[srcPatchi].size(), -1)
                    )
                )
            );
        }
    }

    // Any unset srcPatchFields become calculated
    forAll(srcPatchFields, srcPatchi)
    {
        if (!srcPatchFields.set(srcPatchi))
        {
            // Note: use factory New method instead of direct generation of
            //       calculated so we keep constraints
            srcPatchFields.set
            (
                srcPatchi,
                fvPatchField<Type>::New
                (
                    fvPatchFieldBase::calculatedType(),
                    srcMesh.boundary()[srcPatchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
    }

    auto tresult =
        tmp<VolumeField<Type>>::New
        (
            srcMesh.newIOobject
            (
                IOobject::scopedName
                (
                    type(),
                    "interpolate(" + field.name() + ")"
                )
            ),
            srcMesh,
            field.dimensions(),
            Field<Type>(srcMesh.nCells(), Zero),
            srcPatchFields
        );

    mapTgtToSrc(field, cop, tresult.ref(), secondOrder);

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapTgtToSrc
(
    const tmp<VolumeField<Type>>& tfield,
    const CombineOp& cop,
    const bool secondOrder
) const
{
    return mapTgtToSrc(tfield(), cop, secondOrder);
}


template<class Type>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapTgtToSrc
(
    const VolumeField<Type>& field,
    const bool secondOrder
) const
{
    return mapTgtToSrc(field, plusEqOp<Type>(), secondOrder);
}


template<class Type>
Foam::tmp<Foam::VolumeField<Type>>
Foam::meshToMesh::mapTgtToSrc
(
    const tmp<VolumeField<Type>>& tfield,
    const bool secondOrder
) const
{
    return mapTgtToSrc(tfield(), plusEqOp<Type>(), secondOrder);
}


// ************************************************************************* //
