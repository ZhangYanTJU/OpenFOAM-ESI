/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2025 OpenCFD Ltd.
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

#include "Probes.H"
#include "IOmanip.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SpanStream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Prober>
void Foam::Probes<Prober>::createProbeFiles(const wordList& fieldNames)
{
    // Open new output streams

    bool needsNewFiles = false;
    for (const word& fieldName : fieldNames)
    {
        if (!probeFilePtrs_.found(fieldName))
        {
            needsNewFiles = true;
            break;
        }
    }

    if (needsNewFiles && Pstream::master())
    {
        DebugInfo
            << "Probing fields: " << fieldNames << nl
            << "Probing locations: " << prober_.probeLocations() << nl
            << endl;

        // Put in undecomposed case
        // (Note: gives problems for distributed data running)

        fileName probeDir
        (
            mesh_.time().globalPath()
          / functionObject::outputPrefix
          / name()/mesh_.regionName()
            // Use startTime as the instance for output files
          / mesh_.time().timeName(mesh_.time().startTime().value())
        );
        probeDir.clean();  // Remove unneeded ".."

        // Create directory if needed
        Foam::mkDir(probeDir);

        for (const word& fieldName : fieldNames)
        {
            if (probeFilePtrs_.found(fieldName))
            {
                // Safety
                continue;
            }

            auto osPtr = autoPtr<OFstream>::New(probeDir/fieldName);
            auto& os = *osPtr;

            if (!os.good())
            {
                FatalErrorInFunction
                    << "Cannot open probe output file: " << os.name() << nl
                    << exit(FatalError);
            }

            probeFilePtrs_.insert(fieldName, osPtr);

            DebugInfo<< "open probe stream: " << os.name() << endl;

            const unsigned int width(IOstream::defaultPrecision() + 7);
            os.setf(std::ios_base::left);

            const pointField& probeLocs = prober_.probeLocations();
            const labelList& processors = prober_.processors();
            const labelList& patchIDList = prober_.patchIDList();
            const pointField& oldPoints = prober_.oldPoints();

            forAll(probeLocs, probei)
            {
                os  << "# Probe " << probei << ' ' << probeLocs[probei];

                if (processors[probei] == -1)
                {
                    os  << "  # Not Found";
                }
                else if (probei < patchIDList.size())
                {
                    const label patchi = patchIDList[probei];
                    if (patchi != -1)
                    {
                        const polyBoundaryMesh& bm = mesh_.boundaryMesh();
                        if
                        (
                            patchi < bm.nNonProcessor()
                         || processors[probei] == Pstream::myProcNo()
                        )
                        {
                            os  << " at patch " << bm[patchi].name();
                        }
                        os  << " with a distance of "
                            << mag(probeLocs[probei]-oldPoints[probei])
                            << " m to the original point "
                            << oldPoints[probei];
                    }
                }

                os  << nl;
            }

            os  << setw(width) << "# Time";

            forAll(probeLocs, probei)
            {
                if (prober_.includeOutOfBounds() || processors[probei] != -1)
                {
                    os  << ' ' << setw(width) << probei;
                }
            }
            os  << endl;
        }
    }
}


template<class Prober>
template<class Type>
void Foam::Probes<Prober>::writeValues
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

        OCharStream buf;

        const bool includeOutOfBounds = prober_.includeOutOfBounds();
        const labelList& procs = prober_.processors();
        forAll(values, probei)
        {
            if (includeOutOfBounds || procs[probei] != -1)
            {
                buf.rewind();
                buf << values[probei];
                os  << ' ' << setw(width) << buf.str().data();
            }
        }
        os  << endl;
    }
}


template<class Prober>
template<class GeoField>
void Foam::Probes<Prober>::performAction
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Prober>
Foam::label Foam::Probes<Prober>::prepare(unsigned request)
{
    // Prefilter on selection
    HashTable<wordHashSet> selected =
    (
        loadFromFiles_
      ? IOobjectList(mesh_, mesh_.time().timeName()).classes(fieldSelection_)
      : mesh_.classes(fieldSelection_)
    );

    // Classify and count fields
    label nFields = 0;
    do
    {
        #undef  doLocalCode
        #define doLocalCode(InputType, Target)                                \
        {                                                                     \
            Target.clear();  /* Remove old values */                          \
            const auto iter = selected.cfind(InputType::typeName);            \
            if (iter.good())                                                  \
            {                                                                 \
                /* Add new (current) values */                                \
                Target.append(iter.val().sortedToc());                        \
                nFields += Target.size();                                     \
            }                                                                 \
        }

        doLocalCode(volScalarField, scalarFields_);
        doLocalCode(volVectorField, vectorFields_);
        doLocalCode(volSphericalTensorField, sphericalTensorFields_);
        doLocalCode(volSymmTensorField, symmTensorFields_);
        doLocalCode(volTensorField, tensorFields_);

        doLocalCode(surfaceScalarField, surfaceScalarFields_);
        doLocalCode(surfaceVectorField, surfaceVectorFields_);
        doLocalCode(surfaceSphericalTensorField, surfaceSphericalTensorFields_);
        doLocalCode(surfaceSymmTensorField, surfaceSymmTensorFields_);
        doLocalCode(surfaceTensorField, surfaceTensorFields_);
        #undef doLocalCode
    }
    while (false);


    // Adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields(2*nFields);
        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        DebugInfo
            << "Probing fields: " << currentFields << nl
            << "Probing locations: " << prober_.probeLocations() << nl
            << endl;

        // Close streams for fields that no longer exist
        forAllIters(probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                DebugInfo<< "close probe stream: " << iter()->name() << endl;

                probeFilePtrs_.remove(iter);
            }
        }

        if ((request & ACTION_WRITE) && !currentFields.empty())
        {
            createProbeFiles(currentFields.sortedToc());
        }
    }

    return nFields;
}


template<class Prober>
template<class GeoField>
Foam::tmp<GeoField>
Foam::Probes<Prober>::getOrLoadField(const word& fieldName) const
{
    tmp<GeoField> tfield;

    if (loadFromFiles_)
    {
        tfield.emplace
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            mesh_
        );
    }
    else
    {
        tfield.cref(mesh_.cfindObject<GeoField>(fieldName));
    }

    return tfield;
}


template<class Prober>
template<class Type>
void Foam::Probes<Prober>::storeResults
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Prober>
Foam::Probes<Prober>::Probes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    prober_(mesh_, dict),
    loadFromFiles_(loadFromFiles),
    onExecute_(false),
    fieldSelection_()
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Prober>
bool Foam::Probes<Prober>::verbose(const bool on) noexcept
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


template<class Prober>
bool Foam::Probes<Prober>::read(const dictionary& dict)
{
    dict.readEntry("fields", fieldSelection_);

    verbose_ = dict.getOrDefault("verbose", false);
    onExecute_ = dict.getOrDefault("sampleOnExecute", false);

    // Close old (unused) streams
    prepare(ACTION_NONE);

    return true;
}


template<class Prober>
bool Foam::Probes<Prober>::performAction(unsigned request)
{
    if (!prober_.empty() && request && prepare(request))
    {
        performAction(scalarFields_, request);
        performAction(vectorFields_, request);
        performAction(sphericalTensorFields_, request);
        performAction(symmTensorFields_, request);
        performAction(tensorFields_, request);

        performAction(surfaceScalarFields_, request);
        performAction(surfaceVectorFields_, request);
        performAction(surfaceSphericalTensorFields_, request);
        performAction(surfaceSymmTensorFields_, request);
        performAction(surfaceTensorFields_, request);
    }
    return true;
}


template<class Prober>
bool Foam::Probes<Prober>::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


template<class Prober>
bool Foam::Probes<Prober>::write()
{
    return performAction(ACTION_ALL);
}


// ************************************************************************* //
