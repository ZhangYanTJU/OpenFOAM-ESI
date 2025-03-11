/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
    Foam::SubStrings

Description
    Sub-ranges of a string with a structure similar to std::match_results,
    but without the underlying regular expression matching.

\*---------------------------------------------------------------------------*/

#ifndef Foam_SubStrings_H
#define Foam_SubStrings_H

#include <regex>    // For std::sub_match
#include <memory>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SubStrings Declaration
\*---------------------------------------------------------------------------*/

template<class StringType>
class SubStrings
:
    public std::vector<std::sub_match<typename StringType::const_iterator>>
{
public:

    // Types

        //- The element type
        using value_type =
            typename std::sub_match<typename StringType::const_iterator>;

        //- The const_iterator for the underlying string type
        using string_iterator =
            typename StringType::const_iterator;


    // Member Functions

        //- The total string length of all sub-elements.
        //  Use size() for the number elements.
        std::string::size_type length() const
        {
            std::string::size_type len = 0;

            for (const auto& elem : *this)
            {
                len += elem.length();
            }

            return len;
        }

        //- Return length of element at pos
        std::string::size_type length(size_t pos) const
        {
            return (*this)[pos].length();
        }

        //- Retrieve element at pos, converted to a string type.
        StringType str(size_t pos) const
        {
            return (*this)[pos].str();
        }

        #if __cplusplus >= 201703L
        //- Return element at pos as a string_view
        std::string_view view(size_t pos) const
        {
            #if __cplusplus >= 202002L
            return std::string_view
            (
                (*this)[pos].first,
                (*this)[pos].second
            );
            #else
            // No string_view construct from iterator pairs before c++20
            const auto& s = (*this)[pos];
            return std::string_view
            (
                std::pointer_traits<const char*>::pointer_to(*(s.first)),
                (s.second - s.first)
            );
            #endif
        }
        #endif  /* __cplusplus >= 201703L */

        //- Append sub-string defined by begin/end iterators
        void append
        (
            const typename StringType::const_iterator& b,
            const typename StringType::const_iterator& e
        )
        {
            value_type range;
            range.first = b;
            range.second = e;
            range.matched = true;

            this->push_back(range);
        }


        //- Reduce size by 1 or more elements. Can be called on an empty list.
        void pop_back(size_t n = 1)
        {
            if (n >= this->size())
            {
                this->clear();
            }
            else if (n > 0)
            {
                this->resize(this->size() - n);
            }
        }

        //- Reduce size by 1 or more elements (from the front).
        //- Can be called on an empty list.
        void pop_front(size_t n = 1)
        {
            if (n >= this->size())
            {
                this->clear();
            }
            else if (n > 0)
            {
                // Overlapping range, avoid std::copy, std::move
                for (size_t src = n, dst = 0; src < this->size(); ++src, ++dst)
                {
                    (*this)[dst] = (*this)[src];
                }
                this->resize(this->size() - n);
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
