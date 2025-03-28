/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::faBoundaryMeshEntries

Description
    Read and store dictionary entries for finite-area boundary patches.
    The object is *never* registered to avoid registry name clashes with
    faBoundaryMesh, which may either already have been registered, or
    which should subsequently be registered.

SourceFiles
    faBoundaryMeshEntries.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faBoundaryMeshEntries_H
#define Foam_faBoundaryMeshEntries_H

#include "regIOobject.H"
#include "PtrList.H"
#include "entry.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class faBoundaryMeshEntries Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMeshEntries
:
    public PtrList<entry>,
    public regIOobject
{
public:

    //- Runtime type information
    TypeNameNoDebug("faBoundaryMesh");


    // Constructors

        //- Read construct from IOobject. Never register!
        explicit faBoundaryMeshEntries(const IOobject& io);


    // Factory Methods

        //- Read and return contents. The IOobject is never registered
        static PtrList<entry> readContents(const IOobject& io);


    // Static Functions

        //- Truncate entries at the first processor patch entry
        static void removeProcPatches(PtrList<entry>& entries);

        //- Write list of entries
        static bool writeEntries(Ostream& os, const UPtrList<entry>& entries);

        //- Return a list of patch types, uses the "patch" entry
        static wordList types(const UPtrList<entry>& entries);


   // Member Functions

        //- Truncate at the first processor patch entry
        void removeProcPatches();


    // Characteristics

        //- Return a list of patch types, uses the "patch" entry
        wordList types() const;


    // Write

        //- Write as a plain list of entries
        void writeEntry(Ostream& os) const;

        //- Write as a primitive entry with given name.
        //- If the keyword is empty, revert to a plain list.
        void writeEntry(const word& keyword, Ostream& os) const;

        //- The writeData member function required by regIOobject
        virtual bool writeData(Ostream& os) const;

        //- Write using stream options, forces UNCOMPRESSED
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc = true
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
