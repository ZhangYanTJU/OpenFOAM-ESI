/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "sampledSets.H"
#include "dictionary.H"
#include "Time.H"
#include "globalIndex.H"
#include "volFields.H"
#include "mapPolyMesh.H"
#include "IOmanip.H"
#include "IOobjectList.H"
#include "IndirectList.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSets, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSets,
        dictionary
    );
}

#include "sampledSetsImpl.C"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::coordSetWriter> Foam::sampledSets::newWriter
(
    word writerType,
    const dictionary& topDict,
    const dictionary& setDict
)
{
    // Per-set adjustment
    setDict.readIfPresent<word>("setFormat", writerType);

    return coordSetWriter::New
    (
        writerType,
        // Top-level/set-specific "formatOptions"
        coordSetWriter::formatOptions(topDict, setDict, writerType)
    );
}


Foam::OFstream* Foam::sampledSets::createProbeFile(const word& fieldName)
{
    // Open new output stream

    OFstream* osptr = probeFilePtrs_.lookup(fieldName, nullptr);

    if (!osptr && Pstream::master())
    {
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

        probeFilePtrs_.insert
        (
            fieldName,
            autoPtr<OFstream>::New(probeDir/fieldName)
        );
        osptr = probeFilePtrs_.lookup(fieldName, nullptr);

        if (osptr)
        {
            auto& os = *osptr;

            DebugInfo<< "open probe stream: " << os.name() << endl;

            const unsigned int width(IOstream::defaultPrecision() + 7);

            label nPoints = 0;
            forAll(*this, seti)
            {
                const coordSet& s = gatheredSets_[seti];

                const pointField& pts = static_cast<const pointField&>(s);

                for (const point& p : pts)
                {
                    os  << "# Probe " << nPoints++ << ' ' << p << nl;
                }
            }

            os  << '#' << setw(IOstream::defaultPrecision() + 6)
                << "Probe";

            for (label probei = 0; probei < nPoints; ++probei)
            {
                os  << ' ' << setw(width) << probei;
            }
            os  << nl;

            os  << '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << endl;
        }
    }

    return osptr;
}


void Foam::sampledSets::gatherAllSets()
{
    // Any writer references will become invalid
    for (auto& writer : writers_)
    {
        writer.expire();
    }

    const PtrList<sampledSet>& localSets = *this;

    gatheredSets_.resize_null(localSets.size());
    gatheredSorting_.resize_nocopy(localSets.size());
    globalIndices_.resize_nocopy(localSets.size());

    forAll(localSets, seti)
    {
        const coordSet& coords = localSets[seti];

        globalIndices_[seti].reset(globalIndex::gatherOnly{}, coords.size());
        gatheredSets_.set(seti, coords.gatherSort(gatheredSorting_[seti]));
    }
}


Foam::IOobjectList Foam::sampledSets::preCheckFields(unsigned request)
{
    wordList allFields;    // Just needed for warnings
    HashTable<wordHashSet> selected;

    IOobjectList objects;

    if (loadFromFiles_)
    {
        // Check files for a particular time
        objects = IOobjectList(mesh_, mesh_.time().timeName());

        allFields = objects.names();
        selected = objects.classes(fieldSelection_);
    }
    else
    {
        // Check currently available fields
        allFields = mesh_.names();
        selected = mesh_.classes(fieldSelection_);
    }

    // Parallel consistency (no-op in serial)
    // Probably not needed...
    /// Pstream::mapCombineReduce(selected, HashSetOps::plusEqOp<word>());


    DynamicList<label> missed(fieldSelection_.size());

    // Detect missing fields
    forAll(fieldSelection_, i)
    {
        if
        (
            fieldSelection_[i].isLiteral()
         && !ListOps::found(allFields, fieldSelection_[i])
        )
        {
            missed.append(i);
        }
    }

    if (missed.size() && (request != ACTION_NONE))
    {
        WarningInFunction
            << nl
            << "Cannot find "
            << (loadFromFiles_ ? "field file" : "registered field")
            << " matching "
            << UIndirectList<wordRe>(fieldSelection_, missed) << endl;
    }


    // The selected field names, ordered by (scalar, vector, ...)
    // with internal sorting

    selectedFieldNames_.clear();

    do
    {
        #undef  doLocalCode
        #define doLocalCode(InputType)                                        \
        {                                                                     \
            const auto iter = selected.find(InputType::typeName);             \
            if (iter.good())                                                  \
            {                                                                 \
                selectedFieldNames_.append(iter.val().sortedToc());           \
            }                                                                 \
        }

        doLocalCode(volScalarField);
        doLocalCode(volVectorField);
        doLocalCode(volSphericalTensorField);
        doLocalCode(volSymmTensorField);
        doLocalCode(volTensorField);
        #undef doLocalCode
    }
    while (false);


    // Now propagate field counts (per surface)
    // - can update writer even when not writing without problem

    const label nFields = selectedFieldNames_.size();

    if (writeAsProbes_)
    {
        // Close streams for fields that no longer exist
        forAllIters(probeFilePtrs_, iter)
        {
            if (!selectedFieldNames_.found(iter.key()))
            {
                DebugInfo
                    << "close probe stream: "
                    << iter()->name() << endl;

                probeFilePtrs_.remove(iter);
            }
        }
    }
    else if ((request & ACTION_WRITE) != 0)
    {
        forAll(writers_, seti)
        {
            coordSetWriter& writer = writers_[seti];

            writer.nFields(nFields);
        }
    }

    return objects;
}


void Foam::sampledSets::initDict(const dictionary& dict, const bool initial)
{
    PtrList<sampledSet>::clear();
    if (initial)
    {
        writers_.clear();
    }

    const entry* eptr = dict.findEntry("sets", keyType::LITERAL);

    if (eptr && eptr->isDict())
    {
        PtrList<sampledSet> sampSets(eptr->dict().size());
        if (initial && !writeAsProbes_)
        {
            writers_.resize(sampSets.size());
        }

        label seti = 0;

        for (const entry& dEntry : eptr->dict())
        {
            if (!dEntry.isDict())
            {
                continue;
            }

            const dictionary& subDict = dEntry.dict();

            autoPtr<sampledSet> sampSet =
                sampledSet::New
                (
                    dEntry.keyword(),
                    mesh_,
                    searchEngine_,
                    subDict
                );

            // if (!sampSet || !sampSet->enabled())
            // {
            //     continue;
            // }

            // Define the set
            sampSets.set(seti, sampSet);

            // Define writer, but do not attached
            if (initial && !writeAsProbes_)
            {
                writers_.set
                (
                    seti,
                    newWriter(writeFormat_, dict_, subDict)
                );

                // Use outputDir/TIME/set-name
                writers_[seti].useTimeDir(true);
                writers_[seti].verbose(verbose_);
            }
            ++seti;
        }

        sampSets.resize(seti);
        if (initial && !writeAsProbes_)
        {
            writers_.resize(seti);
        }
        static_cast<PtrList<sampledSet>&>(*this).transfer(sampSets);
    }
    else if (eptr)
    {
        // This is slightly trickier.
        // We want access to the individual dictionaries used for construction

        DynamicList<dictionary> capture;

        PtrList<sampledSet> input
        (
            eptr->stream(),
            sampledSet::iNewCapture(mesh_, searchEngine_, capture)
        );

        PtrList<sampledSet> sampSets(input.size());
        if (initial && !writeAsProbes_)
        {
            writers_.resize(sampSets.size());
        }

        label seti = 0;

        forAll(input, inputi)
        {
            const dictionary& subDict = capture[inputi];

            autoPtr<sampledSet> sampSet = input.release(inputi);

            // if (!sampSet || !sampSet->enabled())
            // {
            //     continue;
            // }

            // Define the set
            sampSets.set(seti, sampSet);

            // Define writer, but do not attached
            if (initial && !writeAsProbes_)
            {
                writers_.set
                (
                    seti,
                    newWriter(writeFormat_, dict_, subDict)
                );

                // Use outputDir/TIME/set-name
                writers_[seti].useTimeDir(true);
                writers_[seti].verbose(verbose_);
            }
            ++seti;
        }

        sampSets.resize(seti);
        if (initial && !writeAsProbes_)
        {
            writers_.resize(seti);
        }

        static_cast<PtrList<sampledSet>&>(*this).transfer(sampSets);
    }

    gatherAllSets();

    needsCorrect_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::sampledSets
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    PtrList<sampledSet>(),
    dict_(dict),
    loadFromFiles_(false),
    verbose_(false),
    onExecute_(false),
    needsCorrect_(false),
    writeAsProbes_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix
      / name/mesh_.regionName()
    ),
    searchEngine_(mesh_),
    samplePointScheme_(),
    writeFormat_(),
    selectedFieldNames_(),
    writers_(),
    probeFilePtrs_(),
    gatheredSets_(),
    gatheredSorting_(),
    globalIndices_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


Foam::sampledSets::sampledSets
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::fvMeshFunctionObject(name, obr, dict),
    PtrList<sampledSet>(),
    dict_(dict),
    loadFromFiles_(loadFromFiles),
    verbose_(false),
    onExecute_(false),
    needsCorrect_(false),
    writeAsProbes_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix
      / name/mesh_.regionName()
    ),
    searchEngine_(mesh_),
    samplePointScheme_(),
    writeFormat_(),
    selectedFieldNames_(),
    writers_(),
    probeFilePtrs_(),
    gatheredSets_(),
    gatheredSorting_(),
    globalIndices_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSets::verbose(const bool on) noexcept
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


bool Foam::sampledSets::read(const dictionary& dict)
{
    if (&dict_ != &dict)
    {
        // Update local copy of dictionary
        dict_ = dict;
    }

    fvMeshFunctionObject::read(dict);

    PtrList<sampledSet>::clear();
    writers_.clear();
    fieldSelection_.clear();
    selectedFieldNames_.clear();

    gatheredSets_.clear();
    gatheredSorting_.clear();
    globalIndices_.clear();

    verbose_ = dict.getOrDefault("verbose", false);
    onExecute_ = dict.getOrDefault("sampleOnExecute", false);

    samplePointScheme_ =
        dict.getOrDefault<word>("interpolationScheme", "cellPoint");

    const entry* eptr = dict.findEntry("sets", keyType::LITERAL);

    if (eptr)
    {
        dict.readEntry("setFormat", writeFormat_);
    }

    // Hard-coded handling of ensemble 'probes' writer
    writeAsProbes_ = ("probes" == writeFormat_);
    if (!writeAsProbes_)
    {
        // Close all streams
        probeFilePtrs_.clear();
    }

    initDict(dict, true);

    // Have some sets, so sort out which fields are needed and report

    if (this->size())
    {
        dict_.readEntry("fields", fieldSelection_);
        fieldSelection_.uniq();

        // Report
        if (writeAsProbes_)
        {
            Info<< "Sampled set as probes ensemble:" << nl;

            forAll(*this, seti)
            {
                const sampledSet& s = (*this)[seti];
                Info<< "  " << s.name();
            }
            Info<< nl;
        }
        else
        {
            Info<< "Sampled set:" << nl;

            forAll(*this, seti)
            {
                const sampledSet& s = (*this)[seti];

                Info<< "    " << s.name() << " -> "
                    << writers_[seti].type() << nl;
            }
        }

        Info<< endl;
    }

    if (debug && Pstream::master())
    {
        Pout<< "sample fields:" << flatOutput(fieldSelection_) << nl
            << "sample sets:" << nl << '(' << nl;

        for
        (
            const sampledSet& s
          : static_cast<const PtrList<sampledSet>&>(*this)
        )
        {
            Pout<< "  " << s << endl;
        }
        Pout<< ')' << endl;
    }

    if (writeAsProbes_)
    {
        (void) preCheckFields(ACTION_NONE);
    }

    // FUTURE:
    // Ensure all sets and merge information are expired
    // expire(true);

    return true;
}


bool Foam::sampledSets::performAction(unsigned request)
{
    if (empty())
    {
        // Nothing to do
        return true;
    }
    else if (needsCorrect_)
    {
        searchEngine_.correct();
        initDict(dict_, false);
    }

    // FUTURE:
    // Update sets and store
    // ...

    // Determine availability of fields.
    // Count number of fields (only seems to be needed for VTK legacy)

    IOobjectList objects = preCheckFields(request);

    const label nFields = selectedFieldNames_.size();

    if (!nFields)
    {
        // Nothing to do
        return true;
    }

    // Update writers
    if (!writeAsProbes_)
    {
        forAll(*this, seti)
        {
            const coordSet& s = gatheredSets_[seti];

            if ((request & ACTION_WRITE) != 0)
            {
                coordSetWriter& writer = writers_[seti];

                if (writer.needsUpdate())
                {
                    writer.setCoordinates(s);
                }

                if (writer.buffering())
                {
                    writer.open
                    (
                        outputPath_
                      / word
                        (
                            s.name()
                          + coordSetWriter::suffix(selectedFieldNames_)
                        )
                    );
                }
                else
                {
                    writer.open(outputPath_/s.name());
                }

                writer.beginTime(mesh_.time());
            }
        }
    }

    // Sample fields

    performAction<VolumeField<scalar>>(objects, request);
    performAction<VolumeField<vector>>(objects, request);
    performAction<VolumeField<sphericalTensor>>(objects, request);
    performAction<VolumeField<symmTensor>>(objects, request);
    performAction<VolumeField<tensor>>(objects, request);


    // Finish this time step
    if (!writeAsProbes_)
    {
        forAll(writers_, seti)
        {
            // Write geometry if no fields were written so that we still
            // can have something to look at

            if ((request & ACTION_WRITE) != 0)
            {
                /// if (!writers_[seti].wroteData())
                /// {
                ///     writers_[seti].write();
                /// }

                writers_[seti].endTime();
            }
        }
    }

    return true;
}


bool Foam::sampledSets::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


bool Foam::sampledSets::write()
{
    return performAction(ACTION_ALL);
}


void Foam::sampledSets::correct()
{
    needsCorrect_ = true;
}


void Foam::sampledSets::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSets::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSets::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //
