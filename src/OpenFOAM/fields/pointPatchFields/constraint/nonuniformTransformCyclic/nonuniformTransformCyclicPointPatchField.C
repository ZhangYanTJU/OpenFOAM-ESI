/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "nonuniformTransformCyclicPointPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonuniformTransformCyclicPointPatchField<Type>::
nonuniformTransformCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(p, iF)
{}


template<class Type>
Foam::nonuniformTransformCyclicPointPatchField<Type>::
nonuniformTransformCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    cyclicPointPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::nonuniformTransformCyclicPointPatchField<Type>::
nonuniformTransformCyclicPointPatchField
(
    const nonuniformTransformCyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    cyclicPointPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::nonuniformTransformCyclicPointPatchField<Type>::
nonuniformTransformCyclicPointPatchField
(
    const nonuniformTransformCyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::nonuniformTransformCyclicPointPatchField<Type>::evaluate
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
        const vectorField& nHat = this->patch().pointNormals();

        const Field<Type> pif(this->patchInternalField());

        // Could write as loop instead...
        tmp<Field<Type>> tvalues
        (
            0.5*(pif + transform(I - 2.0*sqr(nHat), pif))
        );

        // Get internal field to insert values into
        auto& iF = const_cast<Field<Type>&>(this->primitiveField());

        this->setInInternalField(iF, tvalues());
    }
}


// ************************************************************************* //
