/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Description
    Macros to ease declaration of run-time selection tables.

    declareRunTimeSelectionTable is used to create a run-time selection table
    for a base-class which holds constructor pointers on the table.

    declareRunTimeNewSelectionTable is used to create a run-time selection
    table for a derived-class which holds "New" pointers on the table.

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef runTimeSelectionTables_H
#define runTimeSelectionTables_H

#include <memory>   // For std::unique_ptr
#include <utility>  // For std::pair

#include "autoPtr.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Common infrastructure

// Not used directly: declare run-time selection (variables)
#define declareRunTimeSelectionTable0(returnType,prefix,argList)               \
                                                                               \
    /* Function pointer type: construct from argList */                        \
    typedef returnType (*prefix##Ptr)argList;                                  \
                                                                               \
    /* Function table type: construct from argList */                          \
    typedef ::Foam::HashTable                                                  \
    <                                                                          \
        prefix##Ptr,                                                           \
        ::Foam::word,                                                          \
        ::Foam::string::hasher                                                 \
    > prefix##TableType;                                                       \
                                                                               \
    /* Lookup aliases for constructors */                                      \
    typedef ::Foam::HashTable                                                  \
    <                                                                          \
        std::pair<::Foam::word,int>,                                           \
        ::Foam::word,                                                          \
        ::Foam::string::hasher                                                 \
    > prefix##CompatTableType;                                                 \
                                                                               \
    /* Table singleton (storage) */                                            \
    static std::unique_ptr<prefix##TableType> prefix##TablePtr_;               \
                                                                               \
    /* Lookup aliases singleton (storage) */                                   \
    static std::unique_ptr<prefix##CompatTableType> prefix##CompatTablePtr_;   \
                                                                               \
    /* Table singleton (access) */                                             \
    static prefix##TableType& prefix##Table();                                 \
                                                                               \
    /* Aliases singleton (access) */                                           \
    static prefix##CompatTableType& prefix##CompatTable();                     \
                                                                               \
    /* Lookup function pointer */                                              \
    static prefix##Ptr prefix##Table(const ::Foam::word& k)


// Not used directly: storage and access for runtime tables
#define defineRunTimeSelectionTable0(baseType,prefix,Tspecialize)              \
                                                                               \
    /* Define table singleton (storage) */                                     \
    Tspecialize                                                                \
    std::unique_ptr<prefix##TableType> prefix##TablePtr_(nullptr);             \
                                                                               \
    /* Define aliases singleton (storage) */                                   \
    Tspecialize                                                                \
    std::unique_ptr<prefix##CompatTableType> prefix##CompatTablePtr_(nullptr); \
                                                                               \
    /* Define table singleton (access) */                                      \
    Tspecialize prefix##TableType& prefix##Table()                             \
    {                                                                          \
        if (!prefix##TablePtr_)                                                \
        {                                                                      \
            prefix##TablePtr_.reset(new prefix##TableType);                    \
        }                                                                      \
        return *(prefix##TablePtr_);                                           \
    }                                                                          \
                                                                               \
    /* Define table singleton (access) for aliases */                          \
    Tspecialize prefix##CompatTableType& prefix##CompatTable()                 \
    {                                                                          \
        if (!prefix##CompatTablePtr_)                                          \
        {                                                                      \
            prefix##CompatTablePtr_.reset(new prefix##CompatTableType);        \
        }                                                                      \
        return *(prefix##CompatTablePtr_);                                     \
    }                                                                          \
                                                                               \
    /* Define lookup function pointer (singleton) */                           \
    Tspecialize prefix##Ptr prefix##Table(const ::Foam::word& k)               \
    {                                                                          \
        if (prefix##TablePtr_)                                                 \
        {                                                                      \
            const auto& tbl = *prefix##TablePtr_;                              \
            auto iter = tbl.cfind(k);                                          \
            if (!iter.found() && prefix##CompatTablePtr_)                      \
            {                                                                  \
                const auto fnd = prefix##CompatTablePtr_->cfind(k);            \
                if (fnd.found())                                               \
                {                                                              \
                    const auto& alt = fnd.val();                               \
                    iter = tbl.cfind(alt.first);                               \
                    if (::Foam::error::warnAboutAge(alt.second))               \
                    {                                                          \
                        std::cerr                                              \
                            << "Using [v" << alt.second << "] '" << k          \
                            << "' instead of '" << alt.first                   \
                            << " in runtime selection table: " << #baseType    \
                            << std::endl;                                      \
                        ::Foam::error::warnAboutAge("lookup", alt.second);     \
                    }                                                          \
                }                                                              \
            }                                                                  \
            if (iter.found())                                                  \
            {                                                                  \
                return iter.val();                                             \
            }                                                                  \
        }                                                                      \
        return nullptr;                                                        \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Declaration Macros
//

//- Declare run-time selection (variables and adder classes)
#define declareRunTimeSelectionTable\
(AutoPtr,baseType,argNames,argList,parList)                                    \
                                                                               \
    declareRunTimeSelectionTable0(                                             \
        AutoPtr<baseType>, argNames##Constructor, argList);                    \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    struct add##argNames##ConstructorToTable                                   \
    {                                                                          \
        static AutoPtr<baseType> New argList                                   \
        {                                                                      \
            return AutoPtr<baseType>(new baseType##Type parList);              \
        }                                                                      \
                                                                               \
        explicit add##argNames##ConstructorToTable                             \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        {                                                                      \
            if (!argNames##ConstructorTable().insert(k, New))                  \
            {                                                                  \
                std::cerr                                                      \
                    << "Duplicate entry " << k << " in runtime table "         \
                    << #baseType << std::endl;                                 \
                ::Foam::error::safePrintStack(std::cerr);                      \
            }                                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    /* Remove only the entry (not the table) upon destruction */               \
    template<class baseType##Type>                                             \
    struct addRemovable##argNames##ConstructorToTable                          \
    {                                                                          \
        const ::Foam::word name;  /* Lookup name for later removal */          \
                                                                               \
        static AutoPtr<baseType> New argList                                   \
        {                                                                      \
            return AutoPtr<baseType>(new baseType##Type parList);              \
        }                                                                      \
                                                                               \
        explicit addRemovable##argNames##ConstructorToTable                    \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        :                                                                      \
            name(k)                                                            \
        {                                                                      \
            argNames##ConstructorTable().set(k, New);                          \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            if (argNames##ConstructorTablePtr_)                                \
            {                                                                  \
                argNames##ConstructorTablePtr_->erase(name);                   \
            }                                                                  \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
                                                                               \
        void operator=                                                         \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
    };                                                                         \
                                                                               \
    /* Class to add constructor from argList to table */                       \
    template<class baseType##Type>                                             \
    struct addAlias##argNames##ConstructorToTable                              \
    {                                                                          \
        explicit addAlias##argNames##ConstructorToTable                        \
        (                                                                      \
            const ::Foam::word& k,                                             \
            const ::Foam::word& alias,                                         \
            const int ver                                                      \
        )                                                                      \
        {                                                                      \
            argNames##ConstructorCompatTable()                                 \
                .set(alias, std::pair<::Foam::word,int>(k,ver));               \
        }                                                                      \
    };


//- Declare a run-time selection for derived classes
#define declareRunTimeNewSelectionTable\
(AutoPtr,baseType,argNames,argList,parList)                                    \
                                                                               \
    declareRunTimeSelectionTable0(                                             \
        AutoPtr<baseType>,argNames##Constructor,argList);                      \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct add##argNames##ConstructorToTable                                   \
    {                                                                          \
        static AutoPtr<baseType> New##baseType argList                         \
        {                                                                      \
            return AutoPtr<baseType>(baseType##Type::New parList.ptr());       \
        }                                                                      \
                                                                               \
        explicit add##argNames##ConstructorToTable                             \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        {                                                                      \
            if (!argNames##ConstructorTable().insert(k, New##baseType))        \
            {                                                                  \
                std::cerr                                                      \
                    << "Duplicate entry " << k << " in runtime table "         \
                    << #baseType << std::endl;                                 \
                ::Foam::error::safePrintStack(std::cerr);                      \
            }                                                                  \
        }                                                                      \
    };                                                                         \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct addRemovable##argNames##ConstructorToTable                          \
    {                                                                          \
        const ::Foam::word name;  /* Retain name for later removal */          \
                                                                               \
        static AutoPtr<baseType> New##baseType argList                         \
        {                                                                      \
            return AutoPtr<baseType>(baseType##Type::New parList.ptr());       \
        }                                                                      \
                                                                               \
        explicit addRemovable##argNames##ConstructorToTable                    \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        :                                                                      \
            name(k)                                                            \
        {                                                                      \
            argNames##ConstructorTable().set(k, New##baseType);                \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            if (argNames##ConstructorTablePtr_)                                \
            {                                                                  \
                argNames##ConstructorTablePtr_->erase(name);                   \
            }                                                                  \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
        void operator=                                                         \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
    };                                                                         \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct addAlias##argNames##ConstructorToTable                              \
    {                                                                          \
        explicit addAlias##argNames##ConstructorToTable                        \
        (                                                                      \
            const ::Foam::word& k,                                             \
            const ::Foam::word& alias,                                         \
            const int ver                                                      \
        )                                                                      \
        {                                                                      \
            argNames##ConstructorCompatTable()                                 \
                .set(alias, std::pair<::Foam::word,int>(k,ver));               \
        }                                                                      \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Definition Macros
//

//- Define run-time selection table
#define defineRunTimeSelectionTable(baseType,argNames)                         \
                                                                               \
    defineRunTimeSelectionTable0(                                              \
    baseType,baseType::argNames##Constructor,)


//- Define run-time selection table for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateRunTimeSelectionTable(baseType,argNames)                 \
                                                                               \
    defineRunTimeSelectionTable0(                                              \
    baseType,baseType::argNames##Constructor,template<>)


//- Define run-time selection table for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTable(baseType,argNames,Targ)           \
                                                                               \
    defineRunTimeSelectionTable0(                                              \
    baseType,baseType<Targ>::argNames##Constructor,template<>)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
