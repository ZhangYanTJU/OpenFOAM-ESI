/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "uniformFixedValuePointPatchField.H"
#include "SubField.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

// Alternative
// {
//     if (const auto* fpp = isA<facePointPatch>(p).patch())
//     {
//         return fpp->patch();
//     }
//     else
//     {
//         return nullptr;
//     }
// }

template<class Type>
const Foam::polyPatch*
Foam::uniformFixedValuePointPatchField<Type>::getPolyPatch(const pointPatch& p)
{
    const polyMesh& mesh = p.boundaryMesh().mesh()();
    return mesh.boundaryMesh().cfindPatch(p.name());
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF)
{}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{
    if (const polyPatch* pp = this->getPolyPatch(this->patch()))
    {
        refValueFunc_ = PatchFunction1<Type>::New
        (
           *pp,
            "uniformValue",
            dict,
            false  // point values (faceValues = false)
        );
    }
    // Fallback
    refPointValueFunc_ = Function1<Type>::New
    (
        "uniformValue",
        dict,
       &this->internalField().db()
    );

    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        // Evaluate to assign a value
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper)
{
    if (const polyPatch* pp = this->getPolyPatch(this->patch()))
    {
        refValueFunc_ = ptf.refValueFunc_.clone(*pp);
    }
    // Fallback
    refPointValueFunc_ = ptf.refPointValueFunc_.clone();

    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& pfld,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(pfld, iF)
{
    if (const polyPatch* pp = this->getPolyPatch(this->patch()))
    {
        refValueFunc_ = pfld.refValueFunc_.clone(*pp);
    }
    // Fallback
    refPointValueFunc_ = pfld.refPointValueFunc_.clone();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& mapper
)
{
    fixedValuePointPatchField<Type>::autoMap(mapper);

    bool canEvaluate(false);

    if (refValueFunc_)
    {
        refValueFunc_().autoMap(mapper);

        // If mapper is not dependent on time we're ok to evaluate
        if (refValueFunc_->constant())
        {
            canEvaluate = true;
        }
    }
    if (refPointValueFunc_)
    {
        // If mapper is not dependent on time we're ok to evaluate
        if (refPointValueFunc_->constant())
        {
            canEvaluate = true;
        }
    }

    if (canEvaluate)
    {
        this->evaluate();
    }
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchField<Type>::rmap(ptf, addr);

    const auto& tiptf = refCast<const uniformFixedValuePointPatchField>(ptf);

    if (refValueFunc_ && tiptf.refValueFunc_)
    {
        refValueFunc_().rmap(tiptf.refValueFunc_(), addr);
    }
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    const scalar t = this->db().time().timeOutputValue();

    if (refValueFunc_)
    {
        valuePointPatchField<Type>::operator=(refValueFunc_->value(t));
    }
    else
    {
        valuePointPatchField<Type>::operator=(refPointValueFunc_->value(t));
    }
    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::
write(Ostream& os) const
{
    // Note: write value
    fixedValuePointPatchField<Type>::write(os);
    if (refValueFunc_)
    {
        refValueFunc_->writeData(os);
    }
    else if (refPointValueFunc_)
    {
        refPointValueFunc_->writeData(os);
    }
}


// ************************************************************************* //
