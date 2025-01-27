/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2025 OpenCFD Ltd.
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
    Foam::HashFunction

Description
    Hash function class.

    Verify that template overloads are properly resolved

Note
    The second template parameter (bool) is used for SFINAE overloading,

\*---------------------------------------------------------------------------*/

#ifndef HashFunction_H
#define HashFunction_H

#include "Hash.H"

#ifdef FULLDEBUG
#define HashTypeInfo(Name)  \
    static constexpr const char* name() noexcept { return Name; } \
    void info() const { std::cerr<< name() << " hashing\n"; }
#else
#define HashTypeInfo(Name)  void info() const {}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HashFunc Declaration
\*---------------------------------------------------------------------------*/

template<class T>
struct HashFun
{
    void info() const
    {
        #ifdef FULLDEBUG
        if constexpr (std::is_base_of_v<std::string, T>)
        {
            std::cerr<< "std::string hashing\n";
        }
        else
        {
            std::cerr<< "default hashing\n";
        }
        #endif
    }

    unsigned operator()(const T& obj, unsigned seed=0) const
    {
        if constexpr (std::is_base_of_v<std::string, T>)
        {
            return Foam::Hasher(obj.data(), obj.size(), seed);
        }
        else
        {
            return Foam::Hasher(&obj, sizeof(obj), seed);
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing for label
template<> struct HashFun<Foam::label> : Hash<label>
{
    HashTypeInfo("label")
};


//- Hashing for pointers, interpret pointer as a integer type
template<> struct HashFun<void*> : Hash<void *>
{
    HashTypeInfo("pointer")
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Various
#include "edge.H"
#include "face.H"
#include "triFace.H"
#include "Pair.H"
#include "Tuple2.H"
#include "DynamicList.H"
#include "FixedList.H"

namespace Foam
{

template<> struct HashFun<edge> : Hash<edge>
{
    HashTypeInfo("edge")
};


template<> struct HashFun<face> : Hash<face>
{
    HashTypeInfo("face")
};


template<> struct HashFun<triFace> : Hash<triFace>
{
    HashTypeInfo("triFace")
};


template<class T>
struct HashFun<Pair<T>> : Hash<Pair<T>>
{
    HashTypeInfo("Pair")
};

template<class T1, class T2>
struct HashFun<Tuple2<T1, T2>> : Hash<Tuple2<T1, T2>>
{
    HashTypeInfo("Tuple2")
};


template<class T>
struct HashFun<List<T>> : Hash<List<T>>
{
    HashTypeInfo("List")
};

template<class T> struct HashFun<UList<T>> : Hash<UList<T>>
{
    HashTypeInfo("UList")
};

template<class T, int SizeMin>
struct HashFun<DynamicList<T, SizeMin>> : Hash<DynamicList<T, SizeMin>>
{
    HashTypeInfo("DynamicList")
};

template<class T, unsigned N>
struct HashFun<FixedList<T, N>> : Hash<FixedList<T, N>>
{
    HashTypeInfo("FixedList")
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
