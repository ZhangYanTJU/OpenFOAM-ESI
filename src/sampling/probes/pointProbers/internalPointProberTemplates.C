/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "internalPointProber.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalPointProber::sample(const VolumeField<Type>& vField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type>> interpPtr
        (
            interpolation<Type>::New(samplePointScheme_, vField)
        );

        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                const vector& position = operator[](probei);

                values[probei] = interpPtr().interpolate
                (
                    position,
                    elementList_[probei],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                values[probei] = vField[elementList_[probei]];
            }
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalPointProber::sample(const SurfaceField<Type>& sField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    forAll(*this, probei)
    {
        if (faceList_[probei] >= 0)
        {
            values[probei] = sField[faceList_[probei]];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalPointProber::sample(const word& fieldName) const
{
    return sample(mesh_.lookupObject<VolumeField<Type>>(fieldName));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalPointProber::sampleSurfaceField(const word& fieldName) const
{
    return sample(mesh_.lookupObject<SurfaceField<Type>>(fieldName));
}


// ************************************************************************* //
