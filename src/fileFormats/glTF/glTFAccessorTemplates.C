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

#include "glTFAccessor.H"

template<class Type>
Foam::string Foam::glTF::accessor::getValueType() const
{
    string valueType;
    if (pTraits<Type>::typeName == pTraits<label>::typeName)
    {
        valueType = "SCALAR";
    }
    else if (pTraits<Type>::typeName == pTraits<scalar>::typeName)
    {
        valueType = "SCALAR";
    }
    else if (pTraits<Type>::typeName == pTraits<vector>::typeName)
    {
        valueType = "VEC3";
    }
    else if (pTraits<Type>::typeName == pTraits<tensor>::typeName)
    {
        valueType = "MAT3";
    }
    else
    {
        FatalErrorInFunction
            << "Unable to process " << pTraits<Type>::typeName << " fields"
            << abort(FatalError);
    }

    return valueType;
}


template<class Type>
void Foam::glTF::accessor::set(const Field<Type>& fld, bool minMax)
{
    count_ = fld.size();

    type_ = getValueType<Type>();

    componentType_ = key(componentTypes::FLOAT);

    minMax_ = minMax;

    if (minMax_)
    {
        Type minValue = min(fld);
        Type maxValue = max(fld);

        {
            OStringStream ss;
            ss << "[ ";
            for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
            {
                if (dir) ss << ", ";
                ss << float(component(minValue, dir));
            }
            ss << " ]";
            min_ = ss.str();
        }
        {
            OStringStream ss;
            ss << "[ ";
            for (direction dir = 0; dir < pTraits<Type>::nComponents; ++dir)
            {
                if (dir) ss << ", ";
                ss << float(component(maxValue, dir));
            }
            ss << " ]";
            max_ = ss.str();
        }
    }
}


// ************************************************************************* //