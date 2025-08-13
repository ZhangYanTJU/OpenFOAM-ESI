/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2025 OpenCFD Ltd.
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

#include "ProbesBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Prober>
template<class GeoField>
Foam::tmp<GeoField>
Foam::ProbesBase<Prober>::getOrLoadField(const word& fieldName) const
{
    tmp<GeoField> tfield;

    if (loadFromFiles_)
    {
        tfield.reset
        (
            new GeoField
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_
            )
        );
    }
    else
    {
        // Slightly paranoid here
        tfield.cref(mesh_.cfindObject<GeoField>(fieldName));
    }

    return tfield;
}


template<class Prober>
template<class Type>
void Foam::ProbesBase<Prober>::storeResults
(
    const word& fieldName,
    const Field<Type>& values
)
{
    const MinMax<Type> limits(values);
    const Type avgVal = average(values);

    this->setResult("average(" + fieldName + ")", avgVal);
    this->setResult("min(" + fieldName + ")", limits.min());
    this->setResult("max(" + fieldName + ")", limits.max());
    this->setResult("size(" + fieldName + ")", values.size());

    if (verbose_)
    {
        Info<< name() << " : " << fieldName << nl
            << "    avg: " << avgVal << nl
            << "    min: " << limits.min() << nl
            << "    max: " << limits.max() << nl << nl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Prober>
template<class Type>
void Foam::ProbesBase<Prober>::writeValues
(
    const word& fieldName,
    const Field<Type>& values,
    const scalar timeValue
)
{
    if (Pstream::master())
    {
        const unsigned int width(IOstream::defaultPrecision() + 7);
        OFstream& os = *probeFilePtrs_[fieldName];

        os  << setw(width) << timeValue;

        const bool includeOutOfBounds = prober_.includeOutOfBounds();
        const labelList& procs = prober_.processors();
        forAll(values, probei)
        {
            if (includeOutOfBounds || procs[probei] != -1)
            {
                os  << ' ' << setw(width) << values[probei];
            }
        }
        os  << endl;
    }
}


template<class Prober>
template<class GeoField>
void Foam::ProbesBase<Prober>::performAction
(
    const fieldGroup<GeoField>& fieldNames,
    unsigned request
)
{
    for (const word& fieldName : fieldNames)
    {
        tmp<GeoField> tfield = getOrLoadField<GeoField>(fieldName);

        if (tfield)
        {
            const auto& field = tfield();
            const scalar timeValue = field.time().timeOutputValue();

            Field<typename GeoField::value_type> values(prober_.sample(field));

            this->storeResults(fieldName, values);
            if (request & ACTION_WRITE)
            {
                this->writeValues(fieldName, values, timeValue);
            }
        }
    }
}


// ************************************************************************* //
