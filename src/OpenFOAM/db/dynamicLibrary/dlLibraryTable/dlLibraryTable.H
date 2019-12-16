/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

Class
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries.

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class codedBase;

/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Private Data

        DynamicList<void*> libPtrs_;

        DynamicList<fileName> libNames_;


protected:

    // Data Types

        //- Global loader/unloader function type
        //  Called with true on load, false on unload.
        typedef void (*loaderType)(bool);


    // Friends

        //- Allow codedBase as a friend for the low-level hooks
        friend class codedBase;


    // Protected Member Functions

        //- Load/unload hook.
        //  \return true if the function was found and executed
        static bool functionHook
        (
            const bool load,    //!< true = on-load, false = on-unload
            void* handle,
            const std::string& funcName,
            const bool verbose,
            const std::string& context
        );


    // Static Member Functions

        //- Execute global funcName(true) on the library (low-level interface).
        //  Usually done as the first step after opening a library.
        //  \return true if the function was found and executed
        static bool loadHook
        (
            void* handle,
            const std::string& funcName,
            const bool verbose = false,
            const std::string& context = ""
        )
        {
            return functionHook(true, handle, funcName, verbose, context);
        }

        //- Execute global funcName(false) on the library (low-level interface).
        //  Usually done as the last step before closing a library.
        //  \return true if the function was found and executed
        static bool unloadHook
        (
            void* handle,
            const std::string& funcName,
            const bool verbose = false,
            const std::string& context = ""
        )
        {
            return functionHook(false, handle, funcName, verbose, context);
        }


        //- Open specified library name and return pointer.
        //  Warning messages, but no additional side-effects.
        void* openLibrary(const fileName& libName, bool verbose);

        //- No copy construct
        dlLibraryTable(const dlLibraryTable&) = delete;

        //- No copy assignment
        void operator=(const dlLibraryTable&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("dlLibraryTable");

    // Constructors

        //- Construct null
        dlLibraryTable() = default;

        //- Move construct
        dlLibraryTable(dlLibraryTable&&) = default;

        //- Open specified libraries. Ignores duplicate names.
        explicit dlLibraryTable
        (
            const UList<fileName>& libNames,
            bool verbose = true
        );

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary. Verbose = true.
        dlLibraryTable(const dictionary& dict, const word& libsEntry);


    //- Destructor. Closes all libraries loaded by the table.
    ~dlLibraryTable();


    // Member Functions

        //- True if no there are no libraries loaded by the table
        bool empty() const;

        //- The number of libraries loaded by the table
        label size() const;

        //- Clearing closes all libraries loaded by the table.
        void clear(bool verbose = true);

        //- Add to the list of names, but do not yet open.
        //  Ignores duplicate names.
        bool append(const fileName& libName);

        //- Add to the list of names, but do not yet open.
        //  Ignores duplicate names.
        label append(const UList<fileName>& libNames);

        //- Open named, but unopened libraries.
        //- These names will normally have been added with the append() method.
        bool open(bool verbose = true);

        //- Open the named library, optionally warn if problems occur
        void* open(const fileName& libName, bool verbose = true);

        //- Open the named libraries, optionally warn if problems occur
        //  Ignores duplicate names.
        bool open(const UList<fileName>& libNames, bool verbose = true);

        //- Close the named library, optionally warn if problems occur
        bool close(const fileName& libName, bool verbose = true);

        //- Find the handle of the named library
        void* findLibrary(const fileName& libName);

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary.
        bool open(const dictionary& dict, const word& libsEntry);

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary and check the additions
        //- to the given constructor table
        template<class TablePtr>
        bool open
        (
            const dictionary& dict,
            const word& libsEntry,
            const TablePtr& tablePtr
        );


    // Member Operators

        //- Move assignment
        dlLibraryTable& operator=(dlLibraryTable&&) = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
