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

#include "transformFvPatchField.H"
#include "IOstreams.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::transformFvPatchField<Type>::transformFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
   fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::transformFvPatchField<Type>::transformFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{}


template<class Type>
Foam::transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf)
{}


template<class Type>
Foam::transformFvPatchField<Type>::transformFvPatchField
(
    const transformFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // OR (!is_rotational_vectorspace_v<Type>)
    if constexpr (std::is_arithmetic_v<Type>)
    {
        return tmp<Field<Type>>::New(this->size(), pTraits<Type>::one);
    }
    else
    {
        return pTraits<Type>::one - snGradTransformDiag();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        *this
      - cmptMultiply
        (
            valueInternalCoeffs(this->patch().weights()),
            this->patchInternalField()
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFvPatchField<Type>::gradientInternalCoeffs() const
{
    // OR (!is_rotational_vectorspace_v<Type>)
    if constexpr (std::is_arithmetic_v<Type>)
    {
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        return -this->patch().deltaCoeffs()*snGradTransformDiag();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        snGrad()
      - cmptMultiply(gradientInternalCoeffs(), this->patchInternalField());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::transformFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    this->evaluate();
}


// ************************************************************************* //
