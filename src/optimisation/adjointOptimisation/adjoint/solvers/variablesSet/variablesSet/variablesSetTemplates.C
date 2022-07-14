/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "localIOdictionary.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
GeometricField<Type, PatchField, GeoMesh>* variablesSet::allocateNamedField
(
    const fvMesh& mesh,
    const IOobject& io,
    const word& solverName
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;

    // Read-in boundary conditions from given IOobject
    localIOdictionary dict
    (
        IOobject
        (
            io.name(),
            io.instance(),
            io.local(),
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        fieldType::typeName
    );
    dictionary& bField(dict.subDict("boundaryField"));

    // Add solverName to all patch entries.
    // Reduntant if not adjoint fields, but overhead should be small
    for (entry& dEntry : bField)
    {
        if (dEntry.isDict())
        {
            dEntry.dict().add<word>("solverName", solverName, true);
        }
    }
    DebugInfo
        << bField << endl;

    return (new fieldType(io, mesh, dict, true));
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool variablesSet::readFieldOK
(
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& fieldPtr,
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;

    word customName = baseName + solverName;
    IOobject headerCustomName
    (
        customName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    IOobject headerBaseName
    (
        baseName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    bool fieldFound(false);

    // Read field with full name (i.e. baseName plus solverName) if present
    if
    (
        headerCustomName.typeHeaderOk<fieldType>(false)
     && useSolverNameForFields
    )
    {
        fieldPtr.reset
        (
            allocateNamedField<Type, PatchField, GeoMesh>
            (
                mesh,
                headerCustomName,
                solverName
            )
        );
        fieldFound = true;
    }
    // else, see whether the base field exists
    else if (headerBaseName.typeHeaderOk<fieldType>(false))
    {
        fieldPtr.reset
        (
            allocateNamedField<Type, PatchField, GeoMesh>
            (
                mesh,
                headerBaseName,
                solverName
            )
        );

        // Rename field if necessary
        if (useSolverNameForFields)
        {
            Info<< "Field " << customName << " not found" << endl;
            Info<< "Reading base field " << baseName << " and renaming ... "
                << endl;
            fieldPtr.ref().rename(customName);
        }
        fieldFound = true;
    }

    return fieldFound;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::swapAndRename
(
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p1,
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p2
)
{
    // Swapping pointers is OK for the mean flow fields known by the
    // variablesSet (and, in essence, by the solver).
    // The problem is that turbulence models know references to U and phi
    // which cannot be swapped.
    /*
    const word name1 = p1().name();
    const word name2 = p2().name();
    p1.swap(p2);

    p2().rename("temp");
    p1().rename(name1);
    p2().rename(name2);
    */

    // Copy back-up fields to original instead. Slower but there seems to be
    // no other way
    GeometricField<Type, PatchField, GeoMesh> temp("temp", p1());
    p1() == p2();
    p2() == temp;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
variablesSet::allocateRenamedField
(
    const refPtr<GeometricField<Type, PatchField, GeoMesh>>& bf
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    autoPtr<fieldType> returnField(nullptr);
    if (bf.valid())
    {
        returnField = allocateRenamedField(bf());
    }
    return returnField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
variablesSet::allocateRenamedField
(
    const autoPtr<GeometricField<Type, PatchField, GeoMesh>>& bf
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    autoPtr<fieldType> returnField(nullptr);
    if (bf.valid())
    {
        returnField = allocateRenamedField(bf());
    }
    return returnField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
variablesSet::allocateRenamedField
(
    const GeometricField<Type, PatchField, GeoMesh>& bf
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    const word timeName = bf.mesh().time().timeName();
    return autoPtr<fieldType>::New(bf.name() + timeName, bf);
}


template<class Type>
void variablesSet::setField
(
    autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    // Try to read in field with custom or base name
    bool fieldFound
    (
        readFieldOK
        (
            fieldPtr,
            mesh,
            baseName,
            solverName,
            useSolverNameForFields
        )
    );

    // No base or custom field found. This is fatal
    if (!fieldFound)
    {
        FatalErrorInFunction
            << "Could not read field with custom ("
            << word(baseName + solverName) << ") "
            << "or base (" << baseName << ") name"
            << exit(FatalError);
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> variablesSet::allocateField
(
    const fvMesh& mesh,
    const word& baseName,
    const word& solverName,
    const bool useSolverNameForFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    autoPtr<VolFieldType> fieldPtr(nullptr);
    setField(fieldPtr, mesh, baseName, solverName, useSolverNameForFields);

    return tmp<VolFieldType>(fieldPtr.ptr());
}


template<class Type>
void variablesSet::renameTurbulenceField
(
    GeometricField<Type, fvPatchField, volMesh>& baseField,
    const word& solverName
)
{
    // typedefs
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Boundary Boundary;

    // Name of custom field, to be potentially read in
    const word baseName = baseField.name();
    const word customName = baseName + solverName;
    const fvMesh& mesh = baseField.mesh();

    // Renaming of the base field
    baseField.rename(customName);

    // Create field with baseName and write it, to enable continuation
    // Note: gives problems for multi-point runs since we end up with
    // multiple db entries with the same name (one from here and one from
    // the solver that will construct a turbulenceModel).
    // Handled through solver.write() for now
    /*
    if (!mesh.foundObject<VolFieldType>(baseName))
    {
        autoPtr<VolFieldType> baseCopy(new VolFieldType(baseField));
        baseCopy().IOobject::writeOpt(baseField.writeOpt());
        baseCopy().rename(baseName);
        regIOobject::store(baseCopy);
    }
    */

    // Check whether a field with the custom name exists, read it in and
    // set supplied base field to that
    IOobject headerCustomName
    (
        customName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false // Do not register
    );

    if (headerCustomName.typeHeaderOk<VolFieldType>(true))
    {
        Info<< "Reading custom turbulence field " << customName
            << " and replacing " << baseName << nl << endl;
        VolFieldType customField(headerCustomName, mesh);

        // Copy internalfield
        baseField.primitiveFieldRef() = customField.primitiveField();

        // We might apply different boundary conditions per operating point
        // We need to read them from the custom files and substitute the ones
        // known by the turbulence model field
        Boundary& baseBoundary = baseField.boundaryFieldRef();
        Boundary& customBoundary = customField.boundaryFieldRef();
        forAll(baseBoundary, patchI)
        {
            baseBoundary.set
            (
                patchI,
                customBoundary[patchI].clone(baseField.ref())
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::nullifyField
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    field == dimensioned<Type>(field.dimensions(), Zero);
    if (field.nOldTimes())
    {
        fieldType& oldTime = field.oldTime();
        variablesSet::nullifyField(oldTime);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::writeField
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    field.write();
    for (label i = 0; i<field.nOldTimes() - 1; ++i)
    {
        writeField(field.oldTime());
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::setInitField
(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fieldInit,
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    // Store fields at the current and (nOldTimes - 1) time-steps
    // since the last oldTime value is not for restarts
    // Avoiding the storage to a new GeometricField with its oldTimes set
    // due to the overwriting of the latter when time advances
    typedef GeometricField<Type, PatchField, GeoMesh> fieldType;
    fieldInit.setSize(max(field.nOldTimes(), 1));
    fieldInit.set(0, new fieldType(field.name() + "Init" + name(0), field));
    GeometricField<Type, PatchField, GeoMesh>* oldField = &field;
    for (label i = 1; i < fieldInit.size(); ++i)
    {
        oldField = &oldField->oldTime();
        fieldInit.set
        (
            i,
            new fieldType(field.name() + "Init" + name(i), *oldField)
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::restoreFieldInitialization
(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& exact,
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    if (!exact.size())
    {
        FatalErrorInFunction
            << "DataList to be used for restoring primal "
            << "fields is not set" << endl
            << exit(FatalError);
    }
    // Restore GeometricField together with its oldTimes
    field == exact[0];
    label nI(1);
    GeometricField<Type, PatchField, GeoMesh>* oldField = &field;
    // Restore old times stored in exact
    while (nI < exact.size())
    {
        oldField = &oldField->oldTime();
        *oldField == exact[nI];
        ++nI;
    }
    // Covers the rare case of not having read the _0 fields (to be stored in
    // exact) but there is the need of re-setting the oldTimes
    while (nI < field.nOldTimes())
    {
        oldField = &oldField->oldTime();
        *oldField == exact[0];
        nI++;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::resetAdjointField
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    dimensioned<Type> zero(field.dimensions(), Zero);
    field == zero;
    // Restore oldTimes to zero too
    GeometricField<Type, PatchField, GeoMesh>* oldField = &field;
    for (label i = 0; i < field.nOldTimes(); ++i)
    {
        oldField = &oldField->oldTime();
        *oldField == zero;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::setMeanField
(
    refPtr<GeometricField<Type, PatchField, GeoMesh>>& meanFieldPtr,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const fvMesh& mesh
)
{
    meanFieldPtr.reset
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                field.name() + "Mean",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            field,
            false
        )
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void variablesSet::setMeanField
(
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& meanFieldPtr,
    const GeometricField<Type, PatchField, GeoMesh>& field,
    const fvMesh& mesh
)
{
    meanFieldPtr.reset
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                field.name()+"Mean",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            field,
            false
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
