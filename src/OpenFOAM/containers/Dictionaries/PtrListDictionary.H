/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::PtrListDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on the memory managed PtrList
    of \<T\> to provide ordered indexing in addition to the dictionary lookup.

SourceFiles
    PtrListDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrListDictionary_H
#define Foam_PtrListDictionary_H

#include "DictionaryBase.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrListDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrListDictionary
:
    public DictionaryBase<PtrList<T>, T>
{
public:

    //- The template instance used for the dictionary content
    typedef DictionaryBase<PtrList<T>, T> dict_type;


    // Constructors

        //- Default construct: empty without allocation (capacity=0).
        PtrListDictionary() = default;

        //- Construct given initial list size
        explicit PtrListDictionary(const label len)
        {
            dict_type(2*len);
            PtrList<T>::resize(len);
        }

        //- Copy construct
        PtrListDictionary(const PtrListDictionary& dict) = default;

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrListDictionary(Istream& is, const INew& inew)
        :
            dict_type(is, inew)
        {}

        //- Construct from Istream
        explicit PtrListDictionary(Istream& is)
        :
            dict_type(is)
        {}


    // Member Functions

        //- Set element to pointer provided and return old element
        autoPtr<T> set(const label i, const word& key, T* ptr)
        {
            autoPtr<T> old(PtrList<T>::set(i, ptr));

            if (!dict_type::addHashEntry(key, ptr))
            {
                FatalErrorInFunction
                    << "Cannot insert '" << key << "' into hash-table"
                    << Foam::abort(FatalError);
            }
            return old;
        }

        //- Set element to autoPtr value provided and return old element
        autoPtr<T> set(const label i, const word& key, autoPtr<T>& ptr)
        {
            return this->set(i, key, ptr.release());
        }

        //- Set element to tmp value provided and return old element
        autoPtr<T> set(const label i, const word& key, tmp<T>& ptr)
        {
            return this->set(i, key, ptr.ptr());  // release or clone
        }


    // Member Operators

        using PtrList<T>::operator[];

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *(dict_type::lookup(key));
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *(dict_type::lookup(key));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
