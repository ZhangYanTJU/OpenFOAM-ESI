/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "cyclicPeriodicAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicPeriodicAMIFvsPatchField<Type>::cyclicPeriodicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicPeriodicAMIFvsPatchField<Type>::cyclicPeriodicAMIFvsPatchField
(
    const cyclicPeriodicAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p))
{
    if (!isA<cyclicPeriodicAMIFvPatch>(this->patch()))
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
Foam::cyclicPeriodicAMIFvsPatchField<Type>::cyclicPeriodicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    cyclicPeriodicAMIPatch_(refCast<const cyclicPeriodicAMIFvPatch>(p, dict))
{
    if (!isA<cyclicPeriodicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index()
            << " not cyclicPeriodicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicPeriodicAMIFvsPatchField<Type>::cyclicPeriodicAMIFvsPatchField
(
    const cyclicPeriodicAMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    cyclicPeriodicAMIPatch_(ptf.cyclicPeriodicAMIPatch_)
{}


template<class Type>
Foam::cyclicPeriodicAMIFvsPatchField<Type>::cyclicPeriodicAMIFvsPatchField
(
    const cyclicPeriodicAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    cyclicPeriodicAMIPatch_(ptf.cyclicPeriodicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicPeriodicAMIFvsPatchField<Type>::coupled() const
{
    const auto& cpp = this->cyclicPeriodicAMIPatch_;

    if
    (
        Pstream::parRun()
     || (
            cpp.size()
         && cpp.cyclicPeriodicAMIPatch().neighbPatch().size()
        )
    )
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
