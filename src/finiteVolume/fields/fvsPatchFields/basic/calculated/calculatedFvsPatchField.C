/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "calculatedFvsPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF, dict, IOobjectOption::MUST_READ)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedFvsPatchField<Type>::calculatedFvsPatchField
(
    const calculatedFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::fvsPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    auto* patchTypeCtor = patchConstructorTable(p.type());

    if (patchTypeCtor)
    {
        return patchTypeCtor
        (
            p,
            DimensionedField<Type, surfaceMesh>::null()
        );
    }
    else
    {
        return tmp<fvsPatchField<Type>>
        (
            new calculatedFvsPatchField<Type>
            (
                p,
                DimensionedField<Type, surfaceMesh>::null()
            )
        );
    }
}


template<class Type>
template<class AnyType>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::fvsPatchField<Type>::NewCalculatedType
(
    const fvsPatchField<AnyType>& pf
)
{
    return NewCalculatedType(pf.patch());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::calculatedFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    fvsPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
