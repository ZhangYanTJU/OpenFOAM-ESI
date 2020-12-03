/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::IntRange

Description
    An interval of (signed) integers defined by a start and a size.

Note
    Only a minimum of IO operators are defined, to avoid incurring
    too many dependencies or cyclic dependencies.

SourceFiles
    IntRangeI.H
    IntRanges.C

\*---------------------------------------------------------------------------*/

#ifndef IntRange_H
#define IntRange_H

#include "labelFwd.H"
#include <iterator>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
template<class T> class List;

/*---------------------------------------------------------------------------*\
                          Class IntRange Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class IntRange
{
    static_assert(std::is_integral<IntType>::value, "Integral required");

    // Private Data

        //- The start of the interval
        IntType start_;

        //- The length of the interval
        IntType size_;


protected:

    // Protected Member Functions

        //- The value 1 before the start of the range.
        //  Corresponds to the value of the rend() iterator
        inline IntType rend_value() const noexcept;

        //- The value 1 beyond the end of the range.
        //  Corresponds to the value of the end() iterator
        inline IntType end_value() const noexcept;


public:

    // STL type definitions

        //- Type of values the range contains
        typedef IntType value_type;

        //- The type that can represent the size of the range
        typedef IntType size_type;

        //- Input iterator with const access
        class const_iterator;

        //- Reverse input iterator with const access
        class const_reverse_iterator;


    // Generated Methods: copy/move construct, copy/move assignment


    // Constructors

        //- Default construct an empty range (0,0)
        inline constexpr IntRange() noexcept;

        //- Construct a range with specified length, starting at zero (0,len)
        inline explicit constexpr IntRange(const IntType len) noexcept;

        //- Construct a range from start/length, no checks
        inline constexpr IntRange
        (
            const IntType beg,
            const IntType len
        ) noexcept;


    // Member Functions

    // Access

        //- The (inclusive) lower value of the range
        inline IntType start() const noexcept;

        //- The size of the range
        inline IntType size() const noexcept;

        //- Non-const access to start of the range
        inline IntType& start() noexcept;

        //- Non-const access to size of the range
        inline IntType& size() noexcept;

        //- True if range is empty (zero-sized)
        inline bool empty() const noexcept;

        //- The (inclusive) lower value of the range. Same as start()
        inline IntType first() const noexcept;

        //- The (inclusive) upper value of the range
        inline IntType last() const noexcept;


    // Edit

        //- Reset to zero start and zero size
        inline void clear() noexcept;

        //- Reset start and length, no checks
        inline void reset(const IntType beg, const IntType len) noexcept;

        //- Set the start position, no checks
        inline void setStart(const IntType i) noexcept;

        //- Change the size, no checks. Identical to resize()
        inline void setSize(const IntType n) noexcept;

        //- Change the size, no checks. Identical to setSize()
        inline void resize(const IntType n) noexcept;

        //- Enforce non-negative size
        inline void clampSize();


    // Search

        //- True if the (global) value is located within the range
        inline bool found(const IntType value) const noexcept;


    // Member Operators

        //- Offset dereference, without bounds checking
        inline constexpr IntType operator[](const IntType i) const noexcept;

        //- True if the global value is located within the range.
        //  Behaviour identical to found() - usable as a predicate
        inline bool operator()(const IntType value) const noexcept;

        //- Increase the size by 1.
        inline IntType operator++() noexcept;
        inline IntType operator++(int) noexcept;

        //- Increase the size by n.
        inline IntType operator+=(const IntType n) noexcept;

        //- Decrease the size by 1, but never below 0.
        inline IntType operator--() noexcept;
        inline IntType operator--(int) noexcept;

        //- Decrease the size by n, but never below 0.
        inline IntType operator-=(const IntType n) noexcept;

        //- True if range is non-empty
        explicit operator bool() const noexcept { return bool(size_); }


    // Bidirectional input iterators (const)

        //- Return const_iterator to a position within the range,
        //- with bounds checking.
        //  \return iterator at the requested position, or end() for
        //      out-of-bounds
        inline const_iterator at(const IntType i) const;

        //- A const_iterator set to the beginning of the range
        inline const_iterator begin() const;

        //- A const_iterator set to the beginning of the range
        inline const_iterator cbegin() const;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator cend() const;

        //- A const_iterator set to 1 beyond the end of the range.
        inline const_iterator end() const;


    // Bidirectional reverse input iterators (const)

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator rbegin() const;

        //- A const_reverse_iterator set to 1 before the end of range
        inline const_reverse_iterator crbegin() const;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator rend() const;

        //- A const_reverse_iterator set to 1 before the begin of range
        inline const_reverse_iterator crend() const;


    // Iterators

        //- Random-access input iterator with const access
        //  \note No operator+(IntType, iterator) since this provokes
        //      misleading resolution errors
        class const_iterator
        {
            //- The global value
            IntType value_;

        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;
            typedef IntType value_type;
            typedef IntType difference_type;
            typedef const IntType* pointer;
            typedef IntType reference;


        // Constructors

            //- Construct with specified value, or default construct
            inline explicit constexpr const_iterator
            (
                const IntType val = 0
            ) noexcept;


        // Member Operators

            //- Return the value
            inline constexpr IntType operator*() const noexcept;

            //- Offset dereference operator
            inline constexpr IntType operator[](const IntType n) const noexcept;

            //- Prefix increment
            inline const_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_iterator operator++(int);

            //- Prefix decrement
            inline const_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_iterator operator--(int);

            //- Arbitrary increment
            inline const_iterator& operator+=(const IntType n) noexcept;

            //- Arbitrary decrement
            inline const_iterator& operator-=(const IntType n) noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator+
            (
                const IntType n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_iterator operator-
            (
                const IntType n
            ) const noexcept;

            //- Difference operator
            inline constexpr IntType operator-
            (
                const const_iterator& iter
            ) const noexcept;


        // Comparison

            //- Test for equality of values
            inline constexpr bool operator==(const const_iterator& iter)
            const noexcept;

            //- Compare less-than
            inline constexpr bool operator<(const const_iterator& iter)
            const noexcept;


        // Derived comparisons

            constexpr bool operator!=(const const_iterator& iter)
            const noexcept
            {
                return !(*this == iter);
            }

            constexpr bool operator<=(const const_iterator& iter)
            const noexcept
            {
                return !(iter < *this);
            }

            constexpr bool operator>(const const_iterator& iter)
            const noexcept
            {
                return (iter < *this);
            }

            constexpr bool operator>=(const const_iterator& iter)
            const noexcept
            {
                return !(*this < iter);
            }
        };


        //- Random-access reverse input iterator with const access
        //  \note No operator+(IntType, iterator) since this provokes
        //      misleading resolution errors
        class const_reverse_iterator
        {
            //- The global value
            IntType value_;

        public:

            // STL definitions (as per std::iterator)
            typedef std::random_access_iterator_tag iterator_category;
            typedef IntType value_type;
            typedef IntType difference_type;
            typedef const IntType* pointer;
            typedef IntType reference;


        // Constructors

            //- Construct with specified value, or default construct
            inline explicit constexpr const_reverse_iterator
            (
                const IntType val = 0
            ) noexcept;


        // Member Operators

            //- Return the value
            inline constexpr IntType operator*() const noexcept;

            //- Offset dereference operator
            inline constexpr IntType operator[](const IntType n) const noexcept;

            //- Prefix increment
            inline const_reverse_iterator& operator++() noexcept;

            //- Postfix increment
            inline const_reverse_iterator operator++(int);

            //- Prefix decrement
            inline const_reverse_iterator& operator--() noexcept;

            //- Postfix decrement
            inline const_reverse_iterator operator--(int);

            //- Arbitrary increment
            inline const_reverse_iterator& operator+=(const IntType n) noexcept;

            //- Arbitrary decrement
            inline const_reverse_iterator& operator-=(const IntType n) noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator+
            (
                const IntType n
            ) const noexcept;

            //- Return iterator with offset
            inline constexpr const_reverse_iterator operator-
            (
                const IntType n
            ) const noexcept;

            //- Difference operator
            inline constexpr IntType operator-
            (
                const const_reverse_iterator& iter
            )const noexcept;


        // Comparison

            //- Test for equality of values
            inline constexpr bool operator==(const const_reverse_iterator& iter)
            const noexcept;

            //- Reverse compare less-than
            inline constexpr bool operator<(const const_reverse_iterator& iter)
            const noexcept;


        // Derived comparisons

            constexpr bool operator!=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(*this == iter);
            }

            constexpr bool operator<=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(iter < *this);
            }

            constexpr bool operator>(const const_reverse_iterator& iter)
            const noexcept
            {
                return (iter < *this);
            }

            constexpr bool operator>=(const const_reverse_iterator& iter)
            const noexcept
            {
                return !(*this < iter);
            }
        };
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Identity function for common integer types

//- Identity map from an int32_t IntRange
List<label> identity(const IntRange<int32_t>& range);

#if defined(WM_LABEL_SIZE) && (WM_LABEL_SIZE >= 64)
//- Identity map from an int64_t IntRange
List<label> identity(const IntRange<int64_t>& range);
#endif


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Input operators for common integer types

//- Read IntRange from Istream as bracketed (start size) tuple, no checks
Istream& operator>>(Istream& os, IntRange<int32_t>& range);

//- Read IntRange from Istream as bracketed (start size) tuple, no checks
Istream& operator>>(Istream& os, IntRange<int64_t>& range);


// Output operators for common integer types

//- Write IntRange to Ostream as bracketed (start size) tuple
Ostream& operator<<(Ostream& os, const IntRange<int32_t>& range);

//- Write IntRange to Ostream as bracketed (start size) tuple
Ostream& operator<<(Ostream& os, const IntRange<int64_t>& range);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Test for equality of begin/size values
template<class IntType>
inline constexpr bool operator==
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return (a.first() == b.first() && a.size() == b.size());
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size.
template<class IntType>
inline constexpr bool operator<
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         && a.size() < b.size()
        )
    );
}


// Derived comparisons

template<class IntType>
inline constexpr bool operator!=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return !(a == b);
}

template<class IntType>
inline constexpr bool operator<=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return !(b < a);
}

template<class IntType>
inline constexpr bool operator>
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept
{
    return (b < a);
}

template<class IntType>
inline constexpr bool operator>=
(
    const IntRange<IntType>& a,
    const IntRange<IntType>& b
) noexcept

{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IntRangeI.H"

#endif

// ************************************************************************* //