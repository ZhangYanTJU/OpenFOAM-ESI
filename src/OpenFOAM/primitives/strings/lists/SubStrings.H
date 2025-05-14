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

#include <regex>    // For std::ssub_match
#include <memory>
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SubStrings Declaration
\*---------------------------------------------------------------------------*/

class SubStrings
:
    public std::vector<std::ssub_match>
{
public:

    // Static Functions

        //- Return match as a string_view
        static std::string_view view(const std::ssub_match& m)
        {
            if (!m.matched) return std::string_view();

            #if __cplusplus >= 202002L
            return std::string_view(m.first, m.second);
            #else
            // No string_view construct from iterator pairs before c++20
            return std::string_view
            (
                std::pointer_traits<const char*>::pointer_to(*(m.first)),
                (m.second - m.first)
            );
            #endif
        }


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

        //- Retrieve element at pos, as a string
        std::string str(size_t pos) const
        {
            return (*this)[pos].str();
        }

        //- Return element at pos as a string_view
        std::string_view view(size_t pos) const
        {
            return view((*this)[pos]);
        }

        //- Append sub-string defined by begin/end iterators
        void append
        (
            std::string::const_iterator b,
            std::string::const_iterator e
        )
        {
            auto& range = this->emplace_back();
            range.first = b;
            range.second = e;
            range.matched = true;
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
