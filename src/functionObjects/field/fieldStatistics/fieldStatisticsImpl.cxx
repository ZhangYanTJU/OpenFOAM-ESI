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
Foam::functionObjects::fieldStatistics::flatten(const GeoField& fld) const
{
    typedef typename GeoField::value_type value_type;
    typedef Field<value_type> FieldType;

    label n(0);

    if (!internal_)
    {
        // Count boundary values
        for (const auto& pfld : fld.boundaryField())
        {
            if (!pfld.coupled())
            {
                n += pfld.size();
            }
        }
    }

    if (!n)
    {
        // No boundary values - quick return
        return tmp<FieldType>(fld.primitiveField());
    }


    // Combined internal + flattened boundary fields
    // - this adds extra storage, but necessary since the visitor pattern
    //   requires a single input

    auto tflatFld = tmp<FieldType>::New(fld.size() + n);
    auto& flatFld = tflatFld.ref();

    // Copy internal values
    flatFld.slice(0, fld.size()) = fld.primitiveField();

    // Copy boundary values
    n = fld.size();
    for (const auto& pfld : fld.boundaryField())
    {
        if (!pfld.coupled())
        {
            flatFld.slice(n, pfld.size()) = pfld;
            n += pfld.size();
        }
    }

    return tflatFld;
}


template<class Type>
bool Foam::functionObjects::fieldStatistics::calcStat(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const auto* fieldp = obr_.cfindObject<VolFieldType>(fieldName);
    if (!fieldp)
    {
        return false;
    }
    const auto& field = *fieldp;

    tmp<Field<Type>> tfullfield = flatten(field);
    const auto& fullfield = tfullfield.cref();

    HashTable<variantOutput> result;
    for (const auto& iter : statistics_.csorted())
    {
        const statistic& stat = iter.val();

        // Assign a new entry, overwriting existing entries
        result.set(stat.name_, stat.calc(fullfield));
    }

    results_.set(fieldName, result);

    if (extrema_)
    {
        extremaResults_.set
        (
            fieldName,
            (mode_ == modeType::MAG)
          ? calcExtremaData(mag(field)())
          : calcExtremaData(field)
        );
    }

    return true;
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMean(const Field<T>& field) const
{
    if (internal_ && (mean_ == meanType::VOLUMETRIC))
    {
        return gWeightedAverage(mesh_.V(), field);
    }

    return gAverage(field);
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMin(const Field<T>& field) const
{
    return gMin(field);
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcMax(const Field<T>& field) const
{
    return gMax(field);
}


template<class T>
T Foam::functionObjects::fieldStatistics::calcVariance
(
    const Field<T>& field
) const
{
    const T avg(calcMean(field));

    T var = Zero;
    for (const auto& elem : field)
    {
        var += (elem - avg);
    }

    label count = field.size();
    Foam::sumReduce(var, count);

    if (count <= 1)
    {
        return Zero;
    }

    return 1.0/(count - 1.0)*var;
}


template<class GeoField>
Foam::Pair<Foam::functionObjects::fieldStatistics::extremaData>
Foam::functionObjects::fieldStatistics::calcExtremaData
(
    const GeoField& field
) const
{
    typedef typename GeoField::value_type value_type;

    const label proci = UPstream::myProcNo();

    List<value_type> minVs(UPstream::nProcs(), pTraits<value_type>::max);
    List<label> minCells(UPstream::nProcs(), Foam::zero{});
    List<point> minCs(UPstream::nProcs(), Foam::zero{});

    List<value_type> maxVs(UPstream::nProcs(), pTraits<value_type>::min);
    List<label> maxCells(UPstream::nProcs(), Foam::zero{});
    List<point> maxCs(UPstream::nProcs(), Foam::zero{});

    // Find extrema within the internal field
    {
        auto minMaxIds = findMinMax(field);

        if (label celli = minMaxIds.first(); celli >= 0)
        {
            minVs[proci] = field[celli];
            minCells[proci] = celli;
            minCs[proci] = mesh_.C()[celli];
        }

        if (label celli = minMaxIds.second(); celli >= 0)
        {
            maxVs[proci] = field[celli];
            maxCells[proci] = celli;
            maxCs[proci] = mesh_.C()[celli];
        }
    }

    if (!internal_)
    {
        // Find extrema within the boundary fields
        const auto& fieldBoundary = field.boundaryField();
        const auto& CfBoundary = mesh_.C().boundaryField();

        forAll(fieldBoundary, patchi)
        {
            const Field<value_type>& fp = fieldBoundary[patchi];

            if (fp.size())
            {
                const auto& Cfp = CfBoundary[patchi];

                const auto& faceCells
                    = fieldBoundary[patchi].patch().faceCells();

                auto minMaxIds = findMinMax(fp);

                if
                (
                    label facei = minMaxIds.first();
                    facei >= 0 && minVs[proci] > fp[facei]
                )
                {
                    minVs[proci] = fp[facei];
                    minCells[proci] = faceCells[facei];
                    minCs[proci] = Cfp[facei];
                }

                if
                (
                    label facei = minMaxIds.second();
                    facei >= 0 && maxVs[proci] < fp[facei]
                )
                {
                    maxVs[proci] = fp[facei];
                    maxCells[proci] = faceCells[facei];
                    maxCs[proci] = Cfp[facei];
                }
            }
        }
    }

    // Collect info from all processors
    Pstream::allGatherList(minVs);
    Pstream::allGatherList(minCells);
    Pstream::allGatherList(minCs);

    Pstream::allGatherList(maxVs);
    Pstream::allGatherList(maxCells);
    Pstream::allGatherList(maxCs);


    Pair<extremaData> results;

    // min
    {
        auto& slot = results.first();
        const label procId = findMin(minVs);
        slot.value_ = minVs[procId];
        slot.procID_ = procId;
        slot.cellID_ = minCells[procId];
        slot.position_ = minCs[procId];
    }

    // max
    {
        auto& slot = results.second();
        const label procId = findMax(maxVs);
        slot.value_ = maxVs[procId];
        slot.procID_ = procId;
        slot.cellID_ = maxCells[procId];
        slot.position_ = maxCs[procId];
    }

    return results;
}


// ************************************************************************* //
