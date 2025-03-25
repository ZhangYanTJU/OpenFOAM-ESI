/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "pTraits.H"
#include "contiguous.H"
#include "boolVector.H"  // A FixedList pretending to be a vector
#include "complex.H"
#include "vector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "uLabel.H"
#include "Switch.H"

#include <type_traits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

//- Test if Type has typeName member
template<class T, class = void>
struct has_typeName : std::false_type {};

//- Test if Type has typeName member
template<class T>
struct has_typeName
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::typeName)>
>
:
    std::true_type
{};


template<class T>
void printTypeName()
{
    // Both float and double have pTraits typeName = "scalar"!
    if constexpr (std::is_same_v<float, std::remove_cv_t<T>>)
    {
        Info<< "<float>";
    }
    else if constexpr (std::is_same_v<double, std::remove_cv_t<T>>)
    {
        Info<< "<double>";
    }
    else if constexpr (has_typeName<T>::value)
    {
        Info<< pTraits<std::remove_cv_t<T>>::typeName;
    }
    else
    {
        Info<< typeid(T).name();
    }
}


template<class T, class = void>
struct has_zero_one : std::false_type {};

template<class T>
struct has_zero_one
<
    T,
    std::void_t
    <
        decltype(pTraits<std::remove_cv_t<T>>::zero),
        decltype(pTraits<std::remove_cv_t<T>>::one)
    >
> : std::true_type {};


template<class T>
void printMinMaxRange()
{
    if constexpr (has_zero_one<T>::value)
    {
        Info<< " zero=" << pTraits<std::remove_cv_t<T>>::zero
            << " one=" << pTraits<std::remove_cv_t<T>>::one;
    }
}


template<class T>
void printTraits()
{
    printTypeName<T>();
    printMinMaxRange<T>();

    Info<< " integral=" << std::is_integral_v<T>
        << " floating=" << std::is_floating_point_v<T>
        << " rank=" << pTraits_rank<T>::value
        << " nComponents=" << pTraits_nComponents<T>::value
        << " vector-space=" << Switch::name(is_vectorspace_v<T>)
        << " rotate=" << Switch::name(is_rotational_vectorspace_v<T>)
        << " is_label=" << Switch::name(is_contiguous_label<T>::value)
        << " is_scalar=" << Switch::name(is_contiguous_scalar<T>::value)
        << " cmptType=" << typeid(typename pTraits_cmptType<T>::type).name()
        << endl;
}


template<class T>
void printTraits(const pTraits<T>& p)
{
    Info<< p.typeName << " == " << p << endl;
}

template<class T>
void printDecltype()
{
    Info<< "cmptType : " << typeid(T).name() << nl;
}


#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic warning "-Wuninitialized"

int main()
{
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<const complex>();     // Uses specialized pTraits_...
    printTraits<floatVector>();
    printTraits<doubleVector>();
    printTraits<sphericalTensor>();
    printTraits<symmTensor>();
    printTraits<tensor>();
    printTraits<boolVector>();  // Uses specialized pTraits_...
    printTraits<word>();
    printTraits<std::string>();

    Info<< nl;

    {
        pTraits<bool> b(true);
        printTraits(b);
    }

    {
        pTraits<label> l(100);
        printTraits(l);
    }

    printTraits(pTraits<scalar>(3.14159));

    Info<< nl;

    label abc;
    Info<< "uninitialized primitive:"<< abc << endl;

    label def = label();
    Info<< "initialized primitive:"<< def << endl;

    Info<< nl << "some interesting label limits:" << nl;
    std::cout<< "sizeof = " << sizeof(label) << nl;
    std::cout<< "min = " << pTraits<label>::min << nl;
    std::cout<< "max = " << pTraits<label>::max << nl;
    std::cout<< "umax = " << pTraits<uLabel>::max << nl;

    std::cout<< "max_2 = " << pTraits<label>::max/2 << " <=> "
        << (1L << (sizeof(label)*8-2)) << nl;

    std::cout<< "max_4 = " << pTraits<label>::max/4 << " <=> "
        << (1L << (sizeof(label)*8-3)) << nl;

    std::cout<< "max_8 = " << pTraits<label>::max/8 << " <=> "
        << (1L << (sizeof(label)*8-4)) << nl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
