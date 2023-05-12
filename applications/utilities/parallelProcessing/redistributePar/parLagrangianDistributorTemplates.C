/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "parLagrangianDistributor.H"
#include "Time.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "cloud.H"
#include "CompactIOField.H"
#include "DynamicList.H"
#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container>
Foam::wordList Foam::parLagrangianDistributor::filterObjects
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.names<Container>()
      : objects.names<Container>(selectedFields)
    );

    // Parallel synchronise - combine names from all processors
    Pstream::combineReduce(fieldNames, ListOps::uniqueEqOp<word>());
    Foam::sort(fieldNames);  // Consistent order

    return fieldNames;
}


template<class Type>
Foam::label Foam::parLagrangianDistributor::distributeFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef IOField<Type> Container;

    const wordList fieldNames
    (
        filterObjects<IOField<Type>>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject srcIOobject
    (
        "none",
        srcMesh_.time().timeName(),
        cloud::prefix/cloudName,
        srcMesh_,
        IOobject::LAZY_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    IOobject tgtIOobject
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloudName,
        tgtMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    //bool reconstruct = false;

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Distributing lagrangian "
                    << Container::typeName << "s\n" << nl;
            }
            Info<< "        " <<  objectName << nl;
        }
        ++nFields;

        srcIOobject.resetHeader(objectName);
        tgtIOobject.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container field(srcIOobject, haveCloud);

        // Distribute
        map.distribute(field);

        const bool writeOnProc = field.size();

        // Write
        Container(tgtIOobject, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(tgtIOobject.objectPath());
        //}
    }

    if (nFields && verbose_) Info<< endl;

    return nFields;
}


template<class Type>
Foam::label Foam::parLagrangianDistributor::distributeFieldFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef CompactIOField<Field<Type>, Type> Container;

    DynamicList<word> fieldNames;

    // CompactIOField Field names
    fieldNames.push_back
    (
        filterObjects<CompactIOField<Field<Type>, Type>>
        (
            objects,
            selectedFields
        )
    );

    // IOField Field names
    fieldNames.push_back
    (
        filterObjects<IOField<Field<Type>>>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject srcIOobject
    (
        "none",
        srcMesh_.time().timeName(),
        cloud::prefix/cloudName,
        srcMesh_,
        IOobject::LAZY_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    IOobject tgtIOobject
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloudName,
        tgtMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    //bool reconstruct = false;

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Distributing lagrangian "
                    << Container::typeName << "s\n" << nl;
            }
            Info<< "        " <<  objectName << nl;
        }
        ++nFields;

        srcIOobject.resetHeader(objectName);
        tgtIOobject.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container field(srcIOobject, haveCloud);

        // Distribute
        map.distribute(field);

        const bool writeOnProc = field.size();

        // Write
        Container(tgtIOobject, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(tgtIOobject.objectPath());
        //}
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


template<class Container>
Foam::label Foam::parLagrangianDistributor::readFields
(
    const passivePositionParticleCloud& cloud,
    const bool haveCloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    const wordList fieldNames
    (
        filterObjects<Container>
        (
            objects,
            selectedFields
        )
    );

    // Read if present
    IOobject readIO
    (
        "none",
        cloud.time().timeName(),
        cloud,
        IOobject::LAZY_READ,
        IOobject::NO_WRITE,
        IOobject::REGISTER
    );

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reading lagrangian "
                    << Container::typeName << "s\n" << nl;
            }
            Info<< "        " <<  objectName << nl;
        }
        ++nFields;

        readIO.resetHeader(objectName);

        // Read if present (ie, haveCloud means readOnProc)
        Container* fieldPtr = new Container(readIO, haveCloud);

        fieldPtr->store();
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


template<class Container>
Foam::label Foam::parLagrangianDistributor::distributeStoredFields
(
    const mapDistributeBase& map,
    passivePositionParticleCloud& cloud
) const
{
    HashTable<Container*> fields
    (
        cloud.lookupClass<Container>()
    );

    // Parallel-consistent names
    const wordList fieldNames(fields.sortedToc());

    IOobject writeIO
    (
        "none",
        tgtMesh_.time().timeName(),
        cloud::prefix/cloud.name(),
        tgtMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    //bool reconstruct = false;

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        Container& field = *(fields[fieldName]);

        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Distributing lagrangian "
                    << Container::typeName << "s\n" << nl;
            }
            Info<< "        " <<  field.name() << nl;
        }
        ++nFields;

        map.distribute(field);

        writeIO.resetHeader(field.name());

        const bool writeOnProc = field.size();

        // Write
        Container(writeIO, std::move(field)).write(writeOnProc);

        //if (!writeOnProc && !reconstruct)
        //{
        //    // When running with -overwrite it should also delete the old
        //    // files. Below works but is not optimal.
        //
        //    Foam::rm(writeIO.objectPath());
        //}
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


// ************************************************************************* //
