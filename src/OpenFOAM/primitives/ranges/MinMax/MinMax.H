/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::MinMax

Description
    A min/max value pair with additional methods.
    In addition to conveniently storing values, it can be used for logic
    operations or to modify data. A few global functions and functors are
    also provided.

    Examples of use.

    Determine min/max limits from a List of values:
    \verbatim
    List<scalar> values = ...;

    // on construction
    MinMax<scalar> range(values);

    range.reset();

    range += val;

    // global minMax() function
    Info<< minMax(values) << nl;
    \endverbatim

    General comparison operations
    \verbatim
    scalar val;
    if (val < range) ... value is below range min
    if (range.contains(val)) ... value within range
    if (range.compare(val) > 0) ... value is above range max
    if (range(val)) ... value within range - as predicate
    \endverbatim

    Since the range has a predicate form, it can be used as a filter method.
    For example,
    \verbatim
    Info<< "values in range: " << subsetList(values, range) << nl;

    boolList mask = ListOps::create<bool>(values, range);
    Info<< "values values in range " << mask << nl;
    \endverbatim

    One advantage offered by MinMax is to clamp or limit values
    to a particular range. For example,
    \verbatim
    scalarMinMax range(lower, upper);

    scalar val;
    val = range.clamp(val)    .. return clamped values

    // vs.
    val = min(max(value, lower), upper)
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef Foam_MinMax_H
#define Foam_MinMax_H

#include "scalar.H"
#include "zero.H"
#include "Pair.H"
#include "Tuple2.H"
#include "VectorSpace.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class MinMax;

// Common min/max types
typedef MinMax<label> labelMinMax;      //!< A label min/max range
typedef MinMax<scalar> scalarMinMax;    //!< A scalar min/max range


/*---------------------------------------------------------------------------*\
                         Class MinMax Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class MinMax
{
    // Private Data

        //- Min value
        T min_;

        //- Max value
        T max_;

public:

    // Typedefs

        //- The value type the MinMax represents
        typedef T value_type;

        //- Component type
        typedef typename pTraits<T>::cmptType cmptType;


    // Constructors

        //- Default construct: an inverted range
        inline MinMax();

        //- Construct from min/max limits
        inline MinMax(const T& minVal, const T& maxVal);

        //- Implicit construct from min/max limits
        inline MinMax(const std::pair<T,T>& range);

        //- Implicit construct from min/max limits
        inline MinMax(const Pair<T>& range);

        //- Implicit construct from min/max limits
        inline MinMax(const Tuple2<T,T>& range);

        //- Construct with a single zero value
        inline explicit MinMax(Foam::zero);

        //- Implicit construct from zero_one as 0-1 range (pTraits zero, one)
        inline MinMax(Foam::zero_one);

        //- Construct with a single initial value
        inline explicit MinMax(const T& val);

        //- Construct from list of values
        inline explicit MinMax(const UList<T>& vals);


    // Static Member Functions

        //- A semi-infinite range from minVal to the type max
        inline static MinMax<T> ge(const T& minVal);

        //- A semi-infinite range from type min to maxVal
        inline static MinMax<T> le(const T& maxVal);

        //- A 0-1 range corresponding to the pTraits zero, one
        inline static MinMax<T> zero_one();


    // Member Functions

    // Access

        //- The min value
        const T& min() const noexcept { return min_; }

        //- The min value
        T& min() noexcept { return min_; }

        //- The max value
        const T& max() const noexcept { return max_; }

        //- The max value
        T& max() noexcept { return max_; }

        //- The min/max average value
        inline T centre() const;

        //- The min to max span. Zero for invalid range.
        inline T span() const;

        //- The magnitude of the min to max span. Zero for invalid range.
        inline scalar mag() const;

        //- Range is empty if it is inverted
        inline bool empty() const;

        //- Range is non-inverted
        inline bool good() const;

        //- Reset to an inverted (invalid) range
        inline void reset();

        //- Reset min/max to be identical to the specified value
        inline void reset(const T& val);

        //- Reset min/max to specified values
        inline void reset(const T& minVal, const T& maxVal);

        //- Same as reset() - reset to an inverted (invalid) range
        void clear() { reset(); }


    // Testing / Query

        //- Test if the ranges intersect (exclusive check)
        inline bool intersects(const MinMax<T>& b) const;

        //- Test if ranges overlap/touch (inclusive check)
        inline bool overlaps(const MinMax<T>& b) const;

        //- Compares the min/max range with the specified value.
        //  \return
        //  -  0: value is within the range, or range is invalid
        //  - -1: range (max) is less than the value
        //  - +1: range (min) is greater than value
        inline int compare(const T& val) const;

        //- True if the value is within the range (inclusive check)
        inline bool contains(const T& val) const;

        //- Return value clamped component-wise.
        //  If the range is invalid, just returns the value.
        inline T clamp(const T& val) const;


    // Manipulate

        //- Extend the range to include the other min/max range
        inline MinMax<T>& add(const MinMax& other);

        //- Include the value into the range
        inline MinMax<T>& add(const T& val);

        //- Include the values into the range
        inline MinMax<T>& add(const UList<T>& vals);

        //- Include two or more values into the range.
        //  All values must be similar types
        template<class... Args>
        inline MinMax<T>& add(const T& val0, const T& val1, Args&&... values);


    // Member Operators

        //- Identical to contains(), for use as a predicate.
        inline bool operator()(const T& val) const;

        //- Restrict min/max range to union with other range
        inline MinMax<T>& operator&=(const MinMax<T>& b);

        //- Extend min/max range to include other range
        //  Can be used in a reduction operation.
        inline MinMax<T>& operator+=(const MinMax<T>& b);

        //- Extend min/max range to include value
        inline MinMax<T>& operator+=(const T& val);

        //- Extend min/max range to include all values
        inline MinMax<T>& operator+=(const UList<T>& vals);

        //- Multiply range by scalar factor
        inline MinMax<T>& operator*=(scalar s);

        //- Divide range by scalar factor
        inline MinMax<T>& operator/=(scalar s);


    // Housekeeping

        //- Range is non-inverted. Same as good (2022-10)
        bool valid() const { return good(); }

        //- Old method name. Same as clamp (2023-01)
        T clip(const T& val) const { return this->clamp(val); }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Declare MinMax as non-contiguous (similar to Tuple2).
//  Output remains separate (even in binary) and, since the defined
//  \c operator+ is somewhat non-standard, also avoid false matching with
//  any MPI intrinsic operation.
template<class T>
struct is_contiguous<MinMax<T>> : std::false_type {};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Min/max range as a string
template<class T>
word name(const MinMax<T>& range)
{
    return '(' + Foam::name(range.min()) + ',' + Foam::name(range.max()) + ')';
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read min/max range as (min max) tuple from Istream
template<class T>
inline Istream& operator>>(Istream& is, MinMax<T>& range)
{
    is.readBegin("min.max");
    is >> range.min() >> range.max();
    is.readEnd("min.max");

    is.check(FUNCTION_NAME);
    return is;
}

//- Write min/max range as (min max) tuple to Ostream
template<class T>
inline Ostream& operator<<(Ostream& os, const MinMax<T>& range)
{
    os  << token::BEGIN_LIST
        << range.min() << token::SPACE << range.max()
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MinMaxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// Global Functions and Operators
#include "MinMaxOps.H"


// ************************************************************************* //
