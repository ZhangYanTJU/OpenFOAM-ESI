/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2025 OpenCFD Ltd.
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

#include "pointPatchField.H"
#include "pointMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchFieldBase(p),
    internalField_(iF)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchFieldBase(p, dict),
    internalField_(iF)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    pointPatchFieldBase(ptf, p),
    internalField_(iF)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf
)
:
    pointPatchFieldBase(ptf),
    internalField_(ptf.internalField_)
{}


template<class Type>
Foam::pointPatchField<Type>::pointPatchField
(
    const pointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchFieldBase(ptf),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::pointPatchField<Type>::write(Ostream& os) const
{
    os.writeEntry("type", type());

    if (!patchType().empty())
    {
        os.writeEntry("patchType", patchType());
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::patchInternalField
(
    const UList<Type1>& internalData,
    const labelUList& addressing,
    UList<Type1>& pfld
) const
{
    if (FOAM_UNLIKELY(internalData.size() != primitiveField().size()))
    {
        FatalErrorInFunction
            << "Internal field size: " << internalData.size()
            << " != mesh size: " << primitiveField().size() << nl
            << abort(FatalError);
    }

    // For v2412 and earlier this was a field:
    //     const label len = this->size();
    //     pfld.resize_nocopy(len);
    //
    // Now uses pre-sized storage

    const label len = pfld.size();

    #ifdef FULLDEBUG
    if (FOAM_UNLIKELY((addressing.size() < len) || (this->size() < len)))
    {
        FatalErrorInFunction
            << "patchField size = " << len
            << " but patch size = " << this->size()
            << " and addressing size = " << addressing.size() << nl
            << abort(FatalError);
    }
    #endif

    for (label i = 0; i < len; ++i)
    {
        pfld[i] = internalData[addressing[i]];
    }
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::pointPatchField<Type>::patchInternalField
(
    const UList<Type1>& internalData,
    const labelUList& addressing
) const
{
    auto tpfld = tmp<Field<Type1>>::New(this->size());
    this->patchInternalField(internalData, addressing, tpfld.ref());
    return tpfld;
}


template<class Type>
template<class Type1>
Foam::tmp<Foam::Field<Type1>>
Foam::pointPatchField<Type>::patchInternalField
(
    const UList<Type1>& internalData
) const
{
    auto tpfld = tmp<Field<Type1>>::New(this->size());
    this->patchInternalField(internalData, patch().meshPoints(), tpfld.ref());
    return tpfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::pointPatchField<Type>::patchInternalField() const
{
    return patchInternalField(primitiveField());
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    if (FOAM_UNLIKELY(iF.size() != primitiveField().size()))
    {
        FatalErrorInFunction
            << "Internal field size: " << iF.size()
            << " != mesh size: " << primitiveField().size() << nl
            << abort(FatalError);
    }

    if (FOAM_UNLIKELY(pF.size() != size()))
    {
        FatalErrorInFunction
            << "Patch field size: " << pF.size()
            << " != patch size: " << size() << nl
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll(mp, pointi)
    {
        iF[mp[pointi]] += pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelUList& points
) const
{
    if (FOAM_UNLIKELY(iF.size() != primitiveField().size()))
    {
        FatalErrorInFunction
            << "Internal field size: " << iF.size()
            << " != mesh size: " << primitiveField().size() << nl
            << abort(FatalError);
    }

    if (FOAM_UNLIKELY(pF.size() != size()))
    {
        FatalErrorInFunction
            << "Patch field size: " << pF.size()
            << " != patch size: " << size() << nl
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll(points, i)
    {
        label pointi = points[i];
        iF[mp[pointi]] += pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF,
    const labelUList& meshPoints
) const
{
    if (FOAM_UNLIKELY(iF.size() != primitiveField().size()))
    {
        FatalErrorInFunction
            << "Internal field size: " << iF.size()
            << " != mesh size: " << primitiveField().size() << nl
            << abort(FatalError);
    }

    if (FOAM_UNLIKELY(pF.size() != meshPoints.size()))
    {
        FatalErrorInFunction
            << "Patch field size: " << pF.size()
            << " != meshPoints size: " << meshPoints.size() << nl
            << abort(FatalError);
    }

    forAll(meshPoints, pointi)
    {
        iF[meshPoints[pointi]] = pF[pointi];
    }
}


template<class Type>
template<class Type1>
void Foam::pointPatchField<Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    setInInternalField(iF, pF, patch().meshPoints());
}


template<class Type>
void Foam::pointPatchField<Type>::updateCoeffs()
{
    pointPatchFieldBase::setUpdated(true);
}


template<class Type>
void Foam::pointPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!updated())
    {
        updateCoeffs();
    }

    pointPatchFieldBase::setUpdated(false);
    pointPatchFieldBase::setManipulated(false);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pointPatchField<Type>& ptf
)
{
    ptf.write(os);

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
