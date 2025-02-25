/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Simple compilation tests and access for UPstream types

\*---------------------------------------------------------------------------*/

#include "pTraits.H"
#include "contiguous.H"
#include "FixedList.H"
#include "boolVector.H"  // A FixedList pretending to be a vector
#include "barycentric.H"
#include "complex.H"
#include "vector.H"
#include "tensor.H"
#include "uLabel.H"
#include "Switch.H"
#include "IOstreams.H"
#include "UPstream.H"

#include <type_traits>

using namespace Foam;

// Just for debugging
const List<std::string> dataType_names
({
    "byte",
    "int32",
    "int64",
    "uint32",
    "uint64",
    "float",
    "double",
    "long_double",

    "float(2)",
    "double(2)",
    "float(3)",
    "double(3)",
    "float(6)",
    "double(6)",
    "float(9)",
    "double(9)"
});

//- Test for pTraits typeName member : default is false
template<class T, class = void>
struct check_has_typeName : std::false_type {};

//- Test for pTraits zero
template<class T>
struct check_has_typeName
<
    T,
    std::void_t<decltype(pTraits<std::remove_cv_t<T>>::typeName)>
>
:
    std::true_type
{};


// Possible future change...
// //- A supported UPstream data type (intrinsic or user-defined)
// template<>
// struct UPstream_base_dataType<complex> : std::true_type
// {
//     static constexpr auto datatype_id = []()
//     {
//         if constexpr (sizeof(complex) == 2*sizeof(float))
//             return UPstream::dataTypes::type_2float;
//         else
//             return UPstream::dataTypes::type_2double;
//     }();
// };


template<class T>
void printTypeName(const bool showSize = false)
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
    else if constexpr (check_has_typeName<T>::value)
    {
        Info<< pTraits<std::remove_cv_t<T>>::typeName;
    }
    else
    {
        Info<< typeid(T).name();
    }
    if (showSize)
    {
        Info<< " (" << sizeof(T) << " bytes)";
    }
}

template<class Type, bool UseTypeName = true>
void printPstreamTraits(const std::string_view name = std::string_view())
{
    Info<< "========" << nl;
    Info<< "type: ";
    if (!name.empty())
    {
        Info<< name << ' ';
    }
    if constexpr (UseTypeName)
    {
        printTypeName<Type>(true);
    }
    else
    {
        Info<< typeid(Type).name();
        Info<< " (" << sizeof(Type) << " bytes)";
    }

    Info<< ", cmpt:";
    printTypeName<typename Foam::pTraits_cmptType<Type>::type>(true);

    Info<< nl
        << "  is_contiguous:"
        << is_contiguous<Type>::value
        << ", is base:"
        << UPstream_base_dataType<Type>::value
        << ", is cmpt:"
        << UPstream_dataType<Type>::value << nl;

    Info<< "is base:"
        << UPstream_base_dataType<Type>::value
        << " (type:" << int(UPstream_base_dataType<Type>::datatype_id)
        << ")  is alias:" << UPstream_alias_dataType<Type>::value
        << " (type:" << int(UPstream_alias_dataType<Type>::datatype_id)
        << ")" << nl;


    {
        int index = int(UPstream_base_dataType<Type>::datatype_id);
        Info<< "datatype: " << index;

        if (index < dataType_names.size())
        {
            Info<< ' ' << dataType_names[index];
        }
        Info<< nl;
    }

    {
        // Use element or component type (or byte-wise) for data type
        using base = typename UPstream_dataType<Type>::base;
        constexpr auto datatype = UPstream_dataType<Type>::datatype_id;

        Info<< "datatype => ";
        printTypeName<base>();
        Info<< " (" << sizeof(Type)/sizeof(base) << " elems)" << nl
            << "datatype: " << static_cast<int>(datatype) << nl;
    }
}


template<class BinaryOp>
void printOpCodeTraits(BinaryOp bop, std::string_view name)
{
    Info<< "op: " << name << ' ';
    if constexpr (UPstream_opType<BinaryOp>::value)
    {
        Info<< "supported";
    }
    else
    {
        Info<< "unknown";
    }
    Info<< ": " << int(UPstream_opType<BinaryOp>::opcode_id) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    printPstreamTraits<bool>();
    printPstreamTraits<label>();

    printPstreamTraits<int>("<int>");
    printPstreamTraits<long>("<long>");
    printPstreamTraits<unsigned>("<unsigned>");
    printPstreamTraits<unsigned int>("<unsigned int>");
    printPstreamTraits<unsigned long>("<long long>");

    printPstreamTraits<const float>();
    printPstreamTraits<floatVector>();

    printPstreamTraits<scalar>();
    printPstreamTraits<double>();
    printPstreamTraits<doubleVector>();

    // Avoid typeName for barycentric. It is declared, but not defined
    printPstreamTraits<barycentric, false>("barycentric");

    printPstreamTraits<complex>();     // Uses specialized pTraits_...

    printPstreamTraits<boolVector>();  // Uses specialized pTraits_...
    printPstreamTraits<floatVector>();
    printPstreamTraits<doubleVector>();
    printPstreamTraits<tensor>();
    printPstreamTraits<word>();
    printPstreamTraits<std::string>();

    // This will not identify properly at the moment...
    printPstreamTraits<FixedList<doubleVector, 4>>();
    printPstreamTraits<labelPair>();

    Info<< nl
        << "========" << nl
        << "Mapping of binary reduction ops" << nl;

    printOpCodeTraits(minOp<scalar>{}, "min");
    printOpCodeTraits(maxOp<vector>{}, "max");
    printOpCodeTraits(sumOp<vector>{}, "sum");
    printOpCodeTraits(plusOp<vector>{}, "plus");
    printOpCodeTraits(multiplyOp<scalar>{}, "multiply");
    printOpCodeTraits(divideOp<vector>{}, "divide");
    printOpCodeTraits(minMagSqrOp<vector>{}, "minMagSqr");

    printOpCodeTraits(bitAndOp<vector>{}, "bitAnd<vector>");
    printOpCodeTraits(bitOrOp<vector>{}, "bitOr<vector>");

    printOpCodeTraits(bitOrOp<float>{}, "bitOr<float>");
    printOpCodeTraits(bitAndOp<unsigned>{}, "bitAnd<unsigned>");
    printOpCodeTraits(bitOrOp<unsigned>{}, "bitOr<unsigned>");

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
