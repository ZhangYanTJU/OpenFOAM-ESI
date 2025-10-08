/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

Application
    setFields

Group
    grpPreProcessingUtilities

Description
    Set values on a selected set of cells/patch-faces via a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "faMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"
#include "areaFields.H"
#include "coupledFvPatch.H"
#include "coupledFaPatch.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Simple tuple of field type and name (read from stream)
class fieldDescription
{
    word type_;
    word name_;

public:

    const word& type() const noexcept { return type_; }
    const word& name() const noexcept { return name_; }

    explicit fieldDescription(Istream& is)
    {
        is >> type_;
        is >> name_;

        // Eg, read as "volScalarFieldValue", but change to "volScalarField"
        if (type_.ends_with("Value"))
        {
            type_.erase(type_.size()-5);
        }
    }
};


// Consume unused field information
template<class Type>
bool consumeUnusedType(const fieldDescription& fieldDesc, Istream& is)
{
    typedef GeometricField<Type, faPatchField, areaMesh> fieldType1;
    typedef GeometricField<Type, fvPatchField, volMesh>  fieldType2;
    //? typedef GeometricField<Type, faePatchField, areaMesh> fieldType3;
    //? typedef GeometricField<Type, fvsPatchField, volMesh>  fieldType4;

    const bool isExpectedType
    (
        fieldDesc.type() == fieldType1::typeName
     || fieldDesc.type() == fieldType2::typeName
    );

    if (isExpectedType)
    {
        (void) pTraits<Type>(is);
    }

    return isExpectedType;
}


// Consume unused field information
static bool consumeUnused(const fieldDescription& fieldDesc, Istream& is)
{
    return
    (
        consumeUnusedType<scalar>(fieldDesc, is)
     || consumeUnusedType<vector>(fieldDesc, is)
     || consumeUnusedType<sphericalTensor>(fieldDesc, is)
     || consumeUnusedType<symmTensor>(fieldDesc, is)
     || consumeUnusedType<tensor>(fieldDesc, is)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Setting volume fields
template<class Type>
bool setCellFieldType
(
    const fieldDescription& fieldDesc,
    const fvMesh& mesh,
    const labelUList& selectedCells,
    Istream& is
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const bool isExpectedType
    (
        fieldDesc.type() == fieldType::typeName
    );

    if (!isExpectedType)
    {
        return false;
    }

    // Get value from stream
    const Type fieldValue = pTraits<Type>(is);


    // Check the current time directory
    IOobject fieldHeader
    (
        fieldDesc.name(),
        mesh.thisDb().time().timeName(),
        mesh.thisDb(),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    bool found = fieldHeader.typeHeaderOk<fieldType>(true);

    if (!found)
    {
        // Fallback to "constant" directory
        fieldHeader.resetHeader();
        fieldHeader.instance() = mesh.thisDb().time().constant();
        found = fieldHeader.typeHeaderOk<fieldType>(true);
    }

    // Field exists
    if (found)
    {
        Info<< "    - set internal values of "
            << fieldHeader.headerClassName()
            << ": " << fieldDesc.name()
            << " = " << fieldValue << endl;

        fieldType field(fieldHeader, mesh, false);

        if (isNull(selectedCells) || (selectedCells.size() == field.size()))
        {
            field.primitiveFieldRef() = fieldValue;
        }
        else
        {
            for (const label celli : selectedCells)
            {
                field[celli] = fieldValue;
            }
        }

        // Make boundary fields consistent - treat like zeroGradient
        for (auto& pfld : field.boundaryFieldRef())
        {
            pfld = pfld.patchInternalField();
        }

        // Handle any e.g. halo-swaps
        field.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << endl;
        }
    }
    else
    {
        Warning
            << "Field " << fieldDesc.name() << " not found" << endl;
    }

    return isExpectedType;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Setting finite-area face fields
template<class Type>
bool setAreaFieldType
(
    const fieldDescription& fieldDesc,
    const faMesh& mesh,
    const labelUList& selectedFaces,
    Istream& is
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> fieldType;

    const bool isExpectedType
    (
        fieldDesc.type() == fieldType::typeName
    );

    if (!isExpectedType)
    {
        return false;
    }

    // Get value from stream
    const Type fieldValue = pTraits<Type>(is);


    // Check the current time directory
    IOobject fieldHeader
    (
        fieldDesc.name(),
        mesh.thisDb().time().timeName(),
        mesh.thisDb(),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    bool found = fieldHeader.typeHeaderOk<fieldType>(true);

    if (!found)
    {
        // Fallback to "constant" directory
        fieldHeader.resetHeader();
        fieldHeader.instance() = mesh.thisDb().time().constant(),
        found = fieldHeader.typeHeaderOk<fieldType>(true);
    }

    // Field exists
    if (found)
    {
        Info<< "    - set internal values of "
            << fieldHeader.headerClassName()
            << ": " << fieldDesc.name()
            << " = " << fieldValue << endl;

        fieldType field(fieldHeader, mesh);

        if (isNull(selectedFaces) || selectedFaces.size() == field.size())
        {
            field.primitiveFieldRef() = fieldValue;
        }
        else
        {
            for (const label facei : selectedFaces)
            {
                field[facei] = fieldValue;
            }
        }

        // Handle any e.g. halo-swaps
        field.boundaryFieldRef().template evaluateCoupled<coupledFaPatch>();

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << endl;
        }
    }
    else
    {
        Warning
            << "Field " << fieldDesc.name() << " not found" << endl;
    }

    return isExpectedType;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Setting volume boundary fields
template<class Type>
bool setFaceFieldType
(
    const fieldDescription& fieldDesc,
    const fvMesh& mesh,
    const labelUList& selectedFaces,
    Istream& is
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const bool isExpectedType
    (
        fieldDesc.type() == fieldType::typeName
    );

    if (!isExpectedType)
    {
        return false;
    }

    // Get value from stream
    const Type fieldValue = pTraits<Type>(is);


    // Check the current time directory
    IOobject fieldHeader
    (
        fieldDesc.name(),
        mesh.thisDb().time().timeName(),
        mesh.thisDb(),
        IOobjectOption::MUST_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::NO_REGISTER
    );

    bool found = fieldHeader.typeHeaderOk<fieldType>(true);

    if (!found)
    {
        // Fallback to "constant" directory
        fieldHeader.resetHeader();
        fieldHeader.instance() = mesh.thisDb().time().constant();
        found = fieldHeader.typeHeaderOk<fieldType>(true);
    }

    // Field exists
    if (found)
    {
        Info<< "    - set boundary values of "
            << fieldHeader.headerClassName()
            << ": " << fieldDesc.name()
            << " = " << fieldValue << endl;

        fieldType field(fieldHeader, mesh);

        // Create flat list of selected faces and their value.
        Field<Type> allBoundaryValues(mesh.nBoundaryFaces());
        for (const auto& pfld : field.boundaryField())
        {
            SubField<Type>
            (
                allBoundaryValues,
                pfld.size(),
                pfld.patch().offset()
            ) = pfld;
        }

        // Override
        unsigned hasWarned = 0;
        labelList nChanged
        (
            returnReduce(field.boundaryField().size(), maxOp<label>()),
            Zero
        );

        for (const label facei : selectedFaces)
        {
            const label bFacei = facei-mesh.nInternalFaces();

            if (bFacei < 0)
            {
                if (!(hasWarned & 1))
                {
                    hasWarned |= 1;
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << ". Suppressing further warnings." << endl;
                }
            }
            else if (bFacei >= mesh.nBoundaryFaces())
            {
                if (!(hasWarned & 2))
                {
                    hasWarned |= 2;
                    WarningInFunction
                        << "Ignoring out-of-range face " << facei
                        << ". Suppressing further warnings." << endl;
                }
            }
            else
            {
                const label patchi = mesh.boundaryMesh().patchID()[bFacei];

                allBoundaryValues[bFacei] = fieldValue;
                ++nChanged[patchi];
            }
        }

        Pstream::listReduce(nChanged, sumOp<label>());

        auto& fieldBf = field.boundaryFieldRef();

        // Reassign.
        forAll(field.boundaryField(), patchi)
        {
            if (nChanged[patchi] > 0)
            {
                Info<< "    On patch "
                    << field.boundaryField()[patchi].patch().name()
                    << " set " << nChanged[patchi] << " values" << endl;

                fieldBf[patchi] == SubField<Type>
                (
                    allBoundaryValues,
                    fieldBf[patchi].size(),
                    fieldBf[patchi].patch().offset()
                );
            }
        }

        // Handle any e.g. halo-swaps
        field.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << endl;
        }
    }
    else
    {
        Warning
            << "Field " << fieldDesc.name() << " not found" << endl;
    }

    return isExpectedType;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dispatcher for setting volume fields
struct setCellField
{
    autoPtr<setCellField> clone() const { return nullptr; }  // placeholder

    static bool apply
    (
        const fieldDescription& fieldDesc,
        const fvMesh& m,
        const labelUList& selectedCells,
        Istream& is
    )
    {
        return
        (
            setCellFieldType<scalar>(fieldDesc, m, selectedCells, is)
         || setCellFieldType<vector>(fieldDesc, m, selectedCells, is)
         || setCellFieldType<sphericalTensor>(fieldDesc, m, selectedCells, is)
         || setCellFieldType<symmTensor>(fieldDesc, m, selectedCells, is)
         || setCellFieldType<tensor>(fieldDesc, m, selectedCells, is)
        );
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelUList& selected_;

    public:

        iNew(const fvMesh& mesh, const labelUList& selectedCells) noexcept
        :
            mesh_(mesh),
            selected_(selectedCells)
        {}

        autoPtr<setCellField> operator()(Istream& is) const
        {
            const fieldDescription fieldDesc(is);

            bool ok = setCellField::apply(fieldDesc, mesh_, selected_, is);

            if (!ok)
            {
                ok = consumeUnused(fieldDesc, is);

                if (ok)
                {
                    // Not meant for us
                    Info<< "Skip " << fieldDesc.type()
                        << " for finite-volume" << nl;
                }
                else
                {
                    WarningInFunction
                        << "Unsupported field type: "
                        << fieldDesc.type() << endl;
                }
            }

            return nullptr;  // Irrelevant return value
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dispatcher for setting volume boundary fields
struct setFaceField
{
    autoPtr<setFaceField> clone() const { return nullptr; }  // placeholder

    static bool apply
    (
        const fieldDescription& fieldDesc,
        const fvMesh& m,
        const labelUList& selectedFaces,
        Istream& is
    )
    {
        return
        (
            setFaceFieldType<scalar>(fieldDesc, m, selectedFaces, is)
         || setFaceFieldType<vector>(fieldDesc, m, selectedFaces, is)
         || setFaceFieldType<sphericalTensor>(fieldDesc, m, selectedFaces, is)
         || setFaceFieldType<symmTensor>(fieldDesc, m, selectedFaces, is)
         || setFaceFieldType<tensor>(fieldDesc, m, selectedFaces, is)
        );
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelUList& selected_;

    public:

        iNew(const fvMesh& mesh, const labelUList& selectedFaces) noexcept
        :
            mesh_(mesh),
            selected_(selectedFaces)
        {}

        autoPtr<setFaceField> operator()(Istream& is) const
        {
            const fieldDescription fieldDesc(is);

            bool ok = setFaceField::apply(fieldDesc, mesh_, selected_, is);

            if (!ok)
            {
                ok = consumeUnused(fieldDesc, is);

                if (ok)
                {
                    // Not meant for us
                    Info<< "Skip " << fieldDesc.type()
                        << " for finite-volume" << nl;
                }
                else
                {
                    WarningInFunction
                        << "Unsupported field type: "
                        << fieldDesc.type() << endl;
                }
            }

            return nullptr;  // Irrelevant return value
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dispatcher for setting area face fields
struct setAreaField
{
    autoPtr<setAreaField> clone() const { return nullptr; }  // placeholder

    static bool apply
    (
        const fieldDescription& fieldDesc,
        const faMesh& m,
        const labelUList& selectedFaces,
        Istream& is
    )
    {
        return
        (
            setAreaFieldType<scalar>(fieldDesc, m, selectedFaces, is)
         || setAreaFieldType<vector>(fieldDesc, m, selectedFaces, is)
         || setAreaFieldType<sphericalTensor>(fieldDesc, m, selectedFaces, is)
         || setAreaFieldType<symmTensor>(fieldDesc, m, selectedFaces, is)
         || setAreaFieldType<tensor>(fieldDesc, m, selectedFaces, is)
        );
    }

    class iNew
    {
        const faMesh& mesh_;
        const labelUList& selected_;

    public:

        iNew(const faMesh& mesh, const labelUList& selectedFaces) noexcept
        :
            mesh_(mesh),
            selected_(selectedFaces)
        {}

        autoPtr<setAreaField> operator()(Istream& is) const
        {
            const fieldDescription fieldDesc(is);

            bool ok = setAreaField::apply(fieldDesc, mesh_, selected_, is);

            if (!ok)
            {
                ok = consumeUnused(fieldDesc, is);

                if (ok)
                {
                    // Not meant for us
                    Info<< "Skip " << fieldDesc.type()
                        << " for finite-volume" << nl;
                }
                else
                {
                    WarningInFunction
                        << "Unsupported field type: "
                        << fieldDesc.type() << endl;
                }
            }

            return nullptr;  // Irrelevant return value
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();           // Never use function objects

    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addNote
    (
        "Set values on a selected set of cells/patch-faces via a dictionary"
    );

    argList::addOption("dict", "file", "Alternative setFieldsDict");

    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress handling of finite-area mesh/fields"
    );

    #include "addRegionOption.H"
    #include "addAllFaRegionOptions.H"

    // -------------------------

    #include "setRootCase.H"
    #include "createTime.H"

    // Set time from specified time options, or no-op
    timeSelector::setTimeIfPresent(runTime, args);

    #include "createNamedMesh.H"

    // Handle -all-area-regions, -area-regions, -area-region
    #include "getAllFaRegionOptions.H"

    if (args.found("no-finite-area"))
    {
        areaRegionNames.clear();  // For consistency
    }

    // ------------------------------------------------------------------------

    PtrList<faMesh> faMeshes;

    // Setup all area meshes on this region
    if (!areaRegionNames.empty())  // ie, !args.found("no-finite-area")
    {
        faMeshes.resize(areaRegionNames.size());

        label nGoodRegions(0);

        for (const word& areaName : areaRegionNames)
        {
            autoPtr<faMesh> faMeshPtr = faMesh::TryNew(areaName, mesh);

            if (faMeshPtr)
            {
                faMeshes.set(nGoodRegions++, std::move(faMeshPtr));
            }
        }

        faMeshes.resize(nGoodRegions);
    }

    if (faMeshes.size() == 1)
    {
        Info<< "Using finite-area mesh";
        if
        (
            const word& name = polyMesh::regionName(faMeshes[0].name());
            !name.empty()
        )
        {
            Info<< " [" << name << "]";
        }
        Info<< nl;
    }
    else if (faMeshes.size() > 1)
    {
        Info<< "Detected finite-area meshes:";
        for (const faMesh& areaMesh : faMeshes)
        {
            Info<< " [" << areaMesh.name() << "]";
        }
        Info<< nl;
    }

    const word dictName("setFieldsDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary setFieldsDict(dictIO);


    // Note.
    // The PtrList + iNew mechanism is used to trigger the callbacks
    // and perform the desired actions. The contents of the PtrList
    // itself are actually irrelevant.


    // Default field values
    if
    (
        const auto* eptr
      = setFieldsDict.findEntry("defaultFieldValues", keyType::LITERAL)
    )
    {
        ITstream& is = eptr->stream();

        Info<< "Setting volume field default values" << endl;

        PtrList<setCellField> defaultFieldValues
        (
            is,
            setCellField::iNew(mesh, labelList::null())
        );

        for (const faMesh& areaMesh : faMeshes)
        {
            Info<< "Setting area field default values";

            if
            (
                const word& name = polyMesh::regionName(areaMesh.name());
                !name.empty()
            )
            {
                Info<< " [" << name << "]";
            }
            Info<< endl;

            is.rewind();

            PtrList<setAreaField> defaultAreaFieldValues
            (
                is,
                setAreaField::iNew(areaMesh, labelList::null())
            );
        }
        Info<< endl;
    }


    Info<< "Setting field region values" << nl << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    for (const entry& region : regions)
    {
        autoPtr<topoSetSource> source =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        if (source().setType() == topoSetSource::CELLSET_SOURCE)
        {
            cellSet subset
            (
                mesh,
                "cellSet",
                mesh.nCells()/10+1  // Reasonable size estimate.
            );

            source->applyToSet(topoSetSource::NEW, subset);

            labelList selectedCells(subset.sortedToc());

            {
                FixedList<label, 2> stats;
                stats[0] = selectedCells.size();
                stats[1] = mesh.nCells();
                reduce(stats, sumOp<label>());

                Info<< "    Selected " << stats[0] << '/' << stats[1]
                    << " cells" << nl;
            }

            ITstream& is = region.dict().lookup("fieldValues");

            PtrList<setCellField> fieldValues
            (
                is,
                setCellField::iNew(mesh, selectedCells)
            );
        }
        else if (source().setType() == topoSetSource::FACESET_SOURCE)
        {
            faceSet subset
            (
                mesh,
                "faceSet",
                mesh.nBoundaryFaces()/10+1
            );

            source->applyToSet(topoSetSource::NEW, subset);

            labelList selectedFaces(subset.sortedToc());

            Info<< "    Selected " << selectedFaces.size()
                << " faces" << nl;

            ITstream& is = region.dict().lookup("fieldValues");

            PtrList<setFaceField> fieldValues
            (
                is,
                setFaceField::iNew(mesh, selectedFaces)
            );

            for (const faMesh& areaMesh : faMeshes)
            {
                const labelUList& faceLabels = areaMesh.faceLabels();

                // Transcribe from mesh faces to finite-area addressing
                labelList areaFaces(faceLabels.size());

                label nUsed = 0;
                forAll(faceLabels, facei)
                {
                    const label meshFacei = faceLabels[facei];

                    if (subset.test(meshFacei))
                    {
                        areaFaces[nUsed] = facei;
                        ++nUsed;
                    }
                }
                areaFaces.resize(nUsed);

                {
                    FixedList<label, 2> stats;
                    stats[0] = areaFaces.size();
                    stats[1] = faceLabels.size();
                    reduce(stats, sumOp<label>());

                    Info<< "    Selected " << stats[0] << '/' << stats[1]
                        << " area faces for " << areaMesh.name() << endl;
                }

                is.rewind();

                PtrList<setAreaField> fieldValues
                (
                    is,
                    setAreaField::iNew(areaMesh, areaFaces)
                );
            }
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
