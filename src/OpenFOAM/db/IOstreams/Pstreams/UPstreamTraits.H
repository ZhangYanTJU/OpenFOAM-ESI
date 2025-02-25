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
    A set of traits associated with UPstream communication

    - UPstream_dataType trait:
      This wrapper is unwinds the type to check against base/alias, but also
      checks if it is a component aggregate of a supported UPstream data type.
      This will be that main entry point for usage.

    - UPstream_opType trait:
      Provides a mapping of OpenFOAM ops to their MPI equivalent.
      The \c opcode_id is the corresponding internal representation.

Note
    Additional helper traits:

    - UPstream_base_dataType trait:
      Tests true/false if the specified data type has an internal
      MPI equivalent.  The \c datatype_id is the corresponding
      internal enumeration. Even if this tests as false, it will
      always return \c type_byte as the fallback for general contiguous data

    - UPstream_alias_dataType trait:
      Provides mapping for <int/long/long long,...> to the fundamental
      32/64 bit integrals, since <int/long/long long,...> may not otherwise
      directly map on all systems.

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstreamTraits_H
#define Foam_UPstreamTraits_H

#include "UPstream.H"
#include <cstdint>
#include <ios>  // For streamsize
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

// Some vector-space types
// -----------------------
//! \cond
template<class T> class Vector;
template<class T> class SymmTensor;
template<class T> class Tensor;
//! \endcond

// -------------------------
// Some binary operators (as per ops.H), but since ListOps.H is included
// by UPstream.H, don't need to forward declare
// -------------------------
// template<class T> struct minOp;
// template<class T> struct maxOp;
// template<class T> struct plusOp;
// template<class T> struct sumOp;
// template<class T> struct multiplyOp;
// template<class T> struct bitAndOp;
// template<class T> struct bitOrOp;
// template<class T> struct bitXorOp;

//! \cond
template<class T> struct UPstream_dataType;
template<class T> struct UPstream_opType;
//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Base traits

//- A supported UPstream (MPI) reduce/window operation type
template<class T>
struct UPstream_opType : std::false_type
{
    static constexpr auto opcode_id = UPstream::opCodes::invalid;
};


//- A supported UPstream data type (intrinsic or user-defined)
template<class T>
struct UPstream_base_dataType : std::false_type
{
    static constexpr auto datatype_id = UPstream::dataTypes::invalid;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Trait specializations (op-codes)

//- Map minOp\<T\> to \c UPstream::opCodes::op_min
template<class T>
struct UPstream_opType<Foam::minOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_min;
};

//- Map maxOp\<T\> to \c UPstream::opCodes::op_max
template<class T>
struct UPstream_opType<Foam::maxOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_max;
};

//- Map sumOp\<T\> to \c UPstream::opCodes::op_sum
template<class T>
struct UPstream_opType<Foam::sumOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_sum;
};

//- Map plusOp\<T\> to \c UPstream::opCodes::op_sum
//- as a recognized alternative to sumOp\<T\>
template<class T>
struct UPstream_opType<Foam::plusOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_sum;
};

//- Map multiplyOp\<T\> to \c UPstream::opCodes::op_prod
template<class T>
struct UPstream_opType<Foam::multiplyOp<T>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_prod;
};

// NOTE (2025-02):
// currently no mappings provided for
// (op_bool_and, op_bool_or, op_bool_xor) until the calling semantics
// have been properly defined


// These are only viable for unsigned integral types,
// probably not for signed integral types.
// Be extra restrictive for now

//- Map bitAndOp\<T\> to \c UPstream::opCodes::op_bit_and
//- (for unsigned integrals)
template<class T>
struct UPstream_opType<Foam::bitAndOp<T>>
:
    // ie, std::unsigned_integral<T> concept
    std::bool_constant<std::is_integral_v<T> && !std::is_signed_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> && !std::is_signed_v<T>)
            return UPstream::opCodes::op_bit_and;
        else
            return UPstream::opCodes::invalid;
    }();
};

//- Map bitOrOp\<T\> to \c UPstream::opCodes::op_bit_or
//- (for unsigned integrals)
template<class T>
struct UPstream_opType<Foam::bitOrOp<T>>
:
    // ie, std::unsigned_integral<T> concept
    std::bool_constant<std::is_integral_v<T> && !std::is_signed_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> && !std::is_signed_v<T>)
            return UPstream::opCodes::op_bit_or;
        else
            return UPstream::opCodes::invalid;
    }();
};

//- Map bitXorOp\<T\> to \c UPstream::opCodes::op_bit_xor
//- (for unsigned integrals)
template<class T>
struct UPstream_opType<Foam::bitXorOp<T>>
:
    // ie, std::unsigned_integral<T> concept
    std::bool_constant<std::is_integral_v<T> && !std::is_signed_v<T>>
{
    static constexpr auto opcode_id = []() constexpr noexcept
    {
        if constexpr (std::is_integral_v<T> && !std::is_signed_v<T>)
            return UPstream::opCodes::op_bit_xor;
        else
            return UPstream::opCodes::invalid;
    }();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Trait specializations (data types)

// Specializations to match elements of UPstream::dataTypes
#undef  defineUPstreamDataTraits
#define defineUPstreamDataTraits(TypeId, Type)                                \
                                                                              \
    /*! \brief Map \c Type to UPstream::dataTypes::TypeId */                  \
    template<> struct UPstream_base_dataType<Type> : std::true_type           \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };                                                                        \
    /*! \brief Map \c const \c Type to \c UPstream::dataTypes::TypeId */      \
    template<> struct UPstream_base_dataType<const Type> : std::true_type     \
    {                                                                         \
        static constexpr auto datatype_id = UPstream::dataTypes::TypeId;      \
    };


// Intrinsic Types [8]:
// Note: uses 'int32_t,int64_t,...' instead of 'int,long,...' to minimize
// the possibility of duplicates types.
// OpenFOAM defines Foam::label as either int32_t,int64_t (not int,long) too.
defineUPstreamDataTraits(type_byte,   char);
defineUPstreamDataTraits(type_byte,   unsigned char);
defineUPstreamDataTraits(type_int32,  int32_t);
defineUPstreamDataTraits(type_int64,  int64_t);
defineUPstreamDataTraits(type_uint32, uint32_t);
defineUPstreamDataTraits(type_uint64, uint64_t);
defineUPstreamDataTraits(type_float,  float);
defineUPstreamDataTraits(type_double, double);
defineUPstreamDataTraits(type_long_double, long double);

// User Types [6]:
defineUPstreamDataTraits(type_3float,  Vector<float>);
defineUPstreamDataTraits(type_3double, Vector<double>);
defineUPstreamDataTraits(type_6float,  SymmTensor<float>);
defineUPstreamDataTraits(type_6double, SymmTensor<double>);
defineUPstreamDataTraits(type_9float,  Tensor<float>);
defineUPstreamDataTraits(type_9double, Tensor<double>);

#undef defineUPstreamDataTraits


// ------------------------------------------------------------------------- //

//- Explicit handling of data type aliases. This is necessary since
//- different systems map things like 'unsigned long' differently but we
//- restrict ourselves to int32/int64 types
template<class T>
struct UPstream_alias_dataType
:
    std::bool_constant
    <
        // Base type (no alias needed)
        UPstream_base_dataType<std::remove_cv_t<T>>::value ||
        (
            // Or some int 32/64 type to re-map
            std::is_integral_v<T>
         && (sizeof(T) == sizeof(int32_t) || sizeof(T) == sizeof(int64_t))
        )
    >
{
    // Is it using the base type? (no alias needed)
    static constexpr bool is_base =
        UPstream_base_dataType<std::remove_cv_t<T>>::value;

    using base = std::conditional_t
    <
        UPstream_base_dataType<std::remove_cv_t<T>>::value,
        std::remove_cv_t<T>,  // <- using base
        std::conditional_t    // <- using alias
        <
            (
                std::is_integral_v<T>
             && (sizeof(T) == sizeof(int32_t) || sizeof(T) == sizeof(int64_t))
            ),
            std::conditional_t
            <
                (sizeof(T) == sizeof(int32_t)),
                std::conditional_t<std::is_signed_v<T>, int32_t, uint32_t>,
                std::conditional_t<std::is_signed_v<T>, int64_t, uint64_t>
            >,
            char  // Fallback is a byte (eg, arbitrary contiguous data)
        >
    >;

    static constexpr auto datatype_id =
        UPstream_base_dataType<base>::datatype_id;
};


// ------------------------------------------------------------------------- //

//- A supported UPstream data type (fundamental or user-defined)
//- or a component aggregate of a supported UPstream data type.
//
// Is true for the following conditions:
// - The \c Type is directly supported
// - The \c cmptType (eg, from VectorSpace) exists and is directly supported
// - Fallback to byte-wise representation (ie, for contiguous)
// .
template<class T>
struct UPstream_dataType
:
    std::bool_constant
    <
        UPstream_alias_dataType<T>::value
     || UPstream_alias_dataType<typename pTraits_cmptType<T>::type>::value
    >
{
    // Is it using the base type? (ie, not using components)
    static constexpr bool is_base = UPstream_alias_dataType<T>::value;

    //- The underlying data type (if supported) or byte
    using base = std::conditional_t
    <
        UPstream_alias_dataType<T>::value,
        typename UPstream_alias_dataType<T>::base,  // <- using base
        typename UPstream_alias_dataType
        <typename pTraits_cmptType<T>::type>::base  // <- using components
    >;

    //- The corresponding UPstream::dataTypes enumeration
    static constexpr auto datatype_id =
        UPstream_base_dataType<base>::datatype_id;

    //- The size in terms of the number of underlying data elements
    static std::streamsize size(std::streamsize count) noexcept
    {
        if constexpr (UPstream_alias_dataType<T>::value)
        {
            // using base: no multiplier
            return count;
        }
        else
        {
            // using components: with multiplier
            return count*(sizeof(T)/sizeof(base));
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
