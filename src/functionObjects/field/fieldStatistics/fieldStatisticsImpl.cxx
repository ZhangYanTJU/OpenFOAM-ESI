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

#include "fieldStatistics.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<Foam::Field<typename GeoField::value_type>>
Foam::functionObjects::fieldStatistics::flatten(const GeoField& fld)
{
    typedef typename GeoField::value_type value_type;
    typedef Field<value_type> FieldType;


    label n = fld.size();

    if (!internal_)
    {
        for (const auto& pfld : fld.boundaryField())
        {
            if (!pfld.coupled())
            {
                n += pfld.size();
            }
        }
    }

    auto tflatFld = tmp<FieldType>::New(n);
    auto& flatFld = tflatFld.ref();

    // Insert internal values
    SubList<value_type>(flatFld, fld.size(), 0) = fld.primitiveField();

    if (!internal_)
    {
        // Insert boundary values
        n = fld.size();
        for (const auto& pfld : fld.boundaryField())
        {
            if (!pfld.coupled())
            {
                SubList<value_type>(flatFld, pfld.size(), n) = pfld;
                n += pfld.size();
            }
        }
    }

    return tflatFld;
}


template<class Type>
bool Foam::functionObjects::fieldStatistics::calcStat(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (!obr_.foundObject<VolFieldType>(fieldName))
    {
        return false;
    }

    const auto* fieldp = obr_.cfindObject<VolFieldType>(fieldName);
    if (!fieldp)
    {
        return false;
    }
    const auto& field = *fieldp;

    tmp<Field<Type>> tfld = flatten(field);
    const auto& fld = tfld.cref();

    HashTable<variantOutput> result;
    for (const auto& iter : statistics_.csorted())
    {
        const statistic& stat = iter.val();

        // Assign a new entry, overwriting existing entries
        result.set(stat.name_, stat.calc(fld));
    }

    results_.set(fieldName, result);

    return true;
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMean(const Field<T>& field)
{
    if (internal_ && (mean_ == VOLUMETRIC))
    {
        const Field<scalar>& V = mesh_.V();
        return (gSum(V*field)/gSum(V));
    }

    return gAverage(field);
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMin(const Field<T>& field)
{
    const label proci = Pstream::myProcNo();

    // Find min internal field value info

    List<T> minVs(Pstream::nProcs(), pTraits<T>::max);

    labelPair minMaxIds = findMinMax(field);

    label minId = minMaxIds.first();
    if (minId != -1)
    {
        minVs[proci] = field[minId];
    }

    // Collect info from all processors and output
    Pstream::allGatherList(minVs);

    minId = findMin(minVs);

    return minVs[minId];  // minValue
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMax(const Field<T>& field)
{
    const label proci = Pstream::myProcNo();

    // Find max internal field value info

    List<T> maxVs(Pstream::nProcs(), pTraits<T>::min);

    labelPair minMaxIds = findMinMax(field);

    label maxId = minMaxIds.second();
    if (maxId != -1)
    {
        maxVs[proci] = field[maxId];
    }

    // Collect info from all processors and output
    Pstream::allGatherList(maxVs);

    maxId = findMax(maxVs);

    return maxVs[maxId];  // maxValue
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcVariance
(
    const Field<T>& field
)
{
    const T average(calcMean(field));

    T var = Zero;
    for (const auto& elem : field)
    {
        var += (elem - average);
    }

    label n = field.size();

    sumReduce(var, n);

    if (n <= 1)
    {
        return T{};
    }

    return 1/(n-1)*var;
}


// ************************************************************************* //
