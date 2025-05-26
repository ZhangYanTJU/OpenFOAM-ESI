/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "symmetryPlanePointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::symmetryPlanePointPatchField<Type>::symmetryPlanePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    basicSymmetryPointPatchField<Type>(p, iF),
    symmetryPlanePatch_(refCast<const symmetryPlanePointPatch>(p))
{}


template<class Type>
Foam::symmetryPlanePointPatchField<Type>::symmetryPlanePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryPointPatchField<Type>(p, iF, dict),
    symmetryPlanePatch_(refCast<const symmetryPlanePointPatch>(p, dict))
{
    if (!isType<symmetryPlanePointPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not symmetry type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::symmetryPlanePointPatchField<Type>::symmetryPlanePointPatchField
(
    const symmetryPlanePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    basicSymmetryPointPatchField<Type>(ptf, p, iF, mapper),
    symmetryPlanePatch_(refCast<const symmetryPlanePointPatch>(p))
{
    if (!isType<symmetryPlanePointPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::symmetryPlanePointPatchField<Type>::symmetryPlanePointPatchField
(
    const symmetryPlanePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    basicSymmetryPointPatchField<Type>(ptf, iF),
    symmetryPlanePatch_(ptf.symmetryPlanePatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::symmetryPlanePointPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type : no-op
    }
    else
    {
        const symmTensor rot(I - 2.0*sqr(symmetryPlanePatch_.n()));

        const Field<Type> pif(this->patchInternalField());

        // Could write as loop instead...
        tmp<Field<Type>> tvalues
        (
            0.5*(pif + transform(rot, pif))
        );

        // Get internal field to insert values into
        auto& iF = const_cast<Field<Type>&>(this->primitiveField());

        this->setInInternalField(iF, tvalues());
    }
}


// ************************************************************************* //
