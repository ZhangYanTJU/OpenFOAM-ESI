/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fixedNormalSlipPointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    slipPointPatchField<Type>(p, iF),
    n_(vector::max)
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    slipPointPatchField<Type>(p, iF, dict),
    n_(dict.get<vector>("n"))
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const fixedNormalSlipPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    slipPointPatchField<Type>(ptf, p, iF, mapper),
    n_(ptf.n_)
{}


template<class Type>
Foam::fixedNormalSlipPointPatchField<Type>::fixedNormalSlipPointPatchField
(
    const fixedNormalSlipPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    slipPointPatchField<Type>(ptf, iF),
    n_(ptf.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedNormalSlipPointPatchField<Type>::evaluate
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
        tmp<Field<Type>> tvalues
        (
            transform(I - n_*n_, this->patchInternalField())
        );

        // Get internal field to insert values into
        auto& iF = const_cast<Field<Type>&>(this->primitiveField());

        this->setInInternalField(iF, tvalues());
    }
}


template<class Type>
void Foam::fixedNormalSlipPointPatchField<Type>::write(Ostream& os) const
{
    slipPointPatchField<Type>::write(os);
    os.writeEntry("n", n_);
}


// ************************************************************************* //
