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

template<class Type>
bool Foam::functionObjects::fieldStatistics::calcStat(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (!obr_.foundObject<VolFieldType>(fieldName))
    {
        return false;
    }

    const VolFieldType& field = lookupObject<VolFieldType>(fieldName);

    HashTable<variantOutput> result;
    for (const auto& iter : statistics_.csorted())
    {
        const statistic& stat = iter.val();

        // Assign a new entry, overwriting existing entries
        result.set(stat.name_, stat.calc(field));
    }

    results_.set(fieldName, result);

    return true;
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMean(const VolumeField<T>& field)
{
    return gAverage(field);
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMin(const VolumeField<T>& field)
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

    // Find min boundary field info
    const auto& fieldBoundary = field.boundaryField();

    forAll(fieldBoundary, patchi)
    {
        const Field<T>& fp = fieldBoundary[patchi];
        if (fp.size())
        {
            minMaxIds = findMinMax(fp);

            minId = minMaxIds.first();
            if (minVs[proci] > fp[minId])
            {
                minVs[proci] = fp[minId];
            }
        }
    }

    // Collect info from all processors and output
    Pstream::allGatherList(minVs);

    minId = findMin(minVs);

    return minVs[minId];  // minValue
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMax(const VolumeField<T>& field)
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


    // Find max boundary field info
    const auto& fieldBoundary = field.boundaryField();

    forAll(fieldBoundary, patchi)
    {
        const Field<T>& fp = fieldBoundary[patchi];
        if (fp.size())
        {
            minMaxIds = findMinMax(fp);

            maxId = minMaxIds.second();
            if (maxVs[proci] < fp[maxId])
            {
                maxVs[proci] = fp[maxId];
            }
        }
    }

    // Collect info from all processors and output
    Pstream::allGatherList(maxVs);

    maxId = findMax(maxVs);

    return maxVs[maxId];  // maxValue
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcVariance
(
    const VolumeField<T>& field
)
{
    label n = field.size();
    reduce(n, sumOp<label>());

    const T average(calcMean(field));

    T var = Zero;
    for (const auto& elem : field)
    {
        var += (elem - average);
    }
    reduce(var, sumOp<T>());

    if (n <= 1)
    {
        return T{};
    }

    return 1/(n-1)*var;
}


// ************************************************************************* //
