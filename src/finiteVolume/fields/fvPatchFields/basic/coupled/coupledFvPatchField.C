/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF, f)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p, dict)),
    fvPatchField<Type>(p, iF, dict, requireValue)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf.patch())),
    fvPatchField<Type>(ptf)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf.patch())),
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs,
    UList<Type>& result
) const
{
    // Get patch neighbour field, store temporarily in result
    this->patchNeighbourField(result);
    const auto& pnf = result;

    // Same as patchInternalField(...), assuming faceCells are an indirection
    // into internal field, but without additional storage...
    const auto& addr = this->patch().faceCells();
    const auto& iF = this->primitiveField();

    // snGrad = deltaCoeffs * (patchNeighbourField - patchInternalField)

    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = deltaCoeffs[i]*(pnf[i] - iF[addr[i]]);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coupledFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    auto tresult = tmp<Field<Type>>::New(this->size());
    this->snGrad(deltaCoeffs, tresult.ref());
    return tresult;
}


template<class Type>
void Foam::coupledFvPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void Foam::coupledFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        lerp
        (
            this->patchNeighbourField(),
            this->patchInternalField(),
            this->patch().weights()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::coupledFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& tweights,
    UList<Type>& result
) const
{
    const auto& w = tweights();

    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = Type(pTraits<Type>::one)*w[i];
    }
    tweights.clear();
}


template<class Type>
void Foam::coupledFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& tweights,
    UList<Type>& result
) const
{
    const auto& w = tweights();

    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = Type(pTraits<Type>::one)*(1.0 - w[i]);
    }
    tweights.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}


template<class Type>
void Foam::coupledFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs,
    UList<Type>& result
) const
{
    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = -Type(pTraits<Type>::one)*deltaCoeffs[i];
    }
}


template<class Type>
void Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs,
    UList<Type>& result
) const
{
    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = Type(pTraits<Type>::one)*deltaCoeffs[i];
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    auto tresult = tmp<Field<Type>>::New(deltaCoeffs.size());
    this->gradientInternalCoeffs(deltaCoeffs, tresult.ref());
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientInternalCoeffs() const
{
    NotImplemented;
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::coupledFvPatchField<Type>::gradientInternalCoeffs
(
    UList<Type>& result
) const
{
    NotImplemented;
    this->gradientInternalCoeffs(this->patch().deltaCoeffs(), result);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -this->gradientInternalCoeffs(deltaCoeffs);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    NotImplemented;
    return -this->gradientInternalCoeffs();
}


template<class Type>
void Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs
(
    UList<Type>& result
) const
{
    NotImplemented;
    this->gradientBoundaryCoeffs(this->patch().deltaCoeffs(), result);
}


template<class Type>
void Foam::coupledFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
