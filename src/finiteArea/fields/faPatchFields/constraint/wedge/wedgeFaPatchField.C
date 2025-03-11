/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "wedgeFaPatch.H"
#include "wedgeFaPatchField.H"
#include "transformField.H"
#include "symmTransform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFaPatch>(this->patch()))
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
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    this->evaluate();
}


template<class Type>
Foam::wedgeFaPatchField<Type>::wedgeFaPatchField
(
    const wedgeFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::wedgeFaPatchField<Type>::snGrad() const
{
    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type: treat like zero-gradient
        return tmp<Field<Type>>::New(this->size(), Foam::zero{});
    }
    else
    {
        const auto& rot = refCast<const wedgeFaPatch>(this->patch()).faceT();

        const Field<Type> pif(this->patchInternalField());

        const auto& dc = this->patch().deltaCoeffs();

        return
        (
            (0.5*dc)
          * (transform(rot, pif) - pif)
        );
    }
}


template<class Type>
void Foam::wedgeFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if constexpr (!is_rotational_vectorspace_v<Type>)
    {
        // Rotational-invariant type: treat like zero-gradient
        this->extrapolateInternal();
    }
    else
    {
        const auto& rot = refCast<const wedgeFaPatch>(this->patch()).edgeT();

        faPatchField<Type>::operator==
        (
            transform(rot, this->patchInternalField())
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::wedgeFaPatchField<Type>::snGradTransformDiag() const
{
    const auto& rot = refCast<const wedgeFaPatch>(this->patch()).faceT();

    const vector diagV = 0.5*(I - rot).diag();

    return tmp<Field<Type>>::New
    (
        this->size(),
        transformMask<Type>
        (
            pow
            (
                diagV,
                pTraits
                <
                    typename powProduct<vector, pTraits<Type>::rank>::type
                >::zero
            )
        )
    );
}


// ************************************************************************* //
