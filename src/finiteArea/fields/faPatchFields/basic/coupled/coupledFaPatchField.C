/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, f)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    lduInterfaceField(refCast<const lduInterface>(p, dict)),
    faPatchField<Type>(p, iF, dict, requireValue)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf)
{}


template<class Type>
Foam::coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledFaPatchField<Type>::snGrad
(
    UList<Type>& result
) const
{
    // Get patch neighbour field, store temporarily in result
    this->patchNeighbourField(result);
    const auto& pnf = result;

    // Same as patchInternalField(...), assuming edgeFaces are an indirection
    // into internal field, but without additional storage...
    const auto& addr = this->patch().edgeFaces();
    const auto& iF = this->primitiveField();

    const auto& deltaCoeffs = this->patch().deltaCoeffs();

    // snGrad = deltaCoeffs * (patchNeighbourField - patchInternalField)

    const label len = result.size();

    for (label i = 0; i < len; ++i)
    {
        result[i] = deltaCoeffs[i]*(pnf[i] - iF[addr[i]]);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coupledFaPatchField<Type>::snGrad() const
{
    auto tresult = tmp<Field<Type>>::New(this->size());
    this->snGrad(static_cast<UList<Type>&>(tresult.ref()));
    return tresult;
}


template<class Type>
void Foam::coupledFaPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void Foam::coupledFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    Field<Type>::operator=
    (
        lerp
        (
            this->patchNeighbourField(),
            this->patchInternalField(),
            this->patch().weights()
        )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -gradientInternalCoeffs();
}


template<class Type>
void Foam::coupledFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
