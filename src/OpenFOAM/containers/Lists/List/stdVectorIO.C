/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "UList.H"
#include "Istream.H"
#include "Ostream.H"
#include "contiguous.H"
#include "token.H"
#include <vector>

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, std::vector<T>& list)
{
    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("Istream >> std::vector<T> : reading first token");

    if (tok.isCompound())
    {
        // No compound handling ...

        list.clear();  // Clear old contents
        FatalIOErrorInFunction(is)
            << "Support for compoundToken - not implemented" << nl
            << exit(FatalIOError);
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..), int{...} or just a plain '0'

        const label len = tok.labelToken();

        // Resize to length required
        list.resize(len);

        if (is.format() == IOstreamOption::BINARY && is_contiguous<T>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<T>
                (
                    is,
                    reinterpret_cast<char*>(list.data()),   // data_bytes()
                    std::streamsize(list.size())*sizeof(T)  // size_bytes()
                );

                is.fatalCheck
                (
                    "Istream >> std::vector<T> : "
                    "reading binary block"
                );
            }
        }
        else if (std::is_same<char, T>::value)
        {
            // Special treatment for char data (binary I/O only)
            const auto oldFmt = is.format(IOstreamOption::BINARY);

            if (len)
            {
                // read(...) includes surrounding start/end delimiters
                is.read
                (
                    reinterpret_cast<char*>(list.data()),   // data_bytes()
                    std::streamsize(list.size())*sizeof(T)  // size_bytes()
                );

                is.fatalCheck
                (
                    "Istream >> std::vector<char> : "
                    "reading binary block"
                );
            }

            is.format(oldFmt);
        }
        else
        {
            // Begin of contents marker
            const char delimiter = is.readBeginList("List");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    auto iter = list.begin();
                    const auto last = list.end();

                    // Contents
                    for (/*nil*/; (iter != last); (void)++iter)
                    {
                        is >> *iter;

                        is.fatalCheck
                        (
                            "Istream >> std::vector<char> : "
                            "reading entry"
                        );
                    }
                }
                else
                {
                    // Uniform content (delimiter == token::BEGIN_BLOCK)

                    T elem;
                    is >> elem;

                    is.fatalCheck
                    (
                        "Istream >> std::vector<char> : "
                        "reading the single entry"
                    );

                    // Fill with the value
                    list.assign(list.size(), elem);
                }
            }

            // End of contents marker
            is.readEndList("List");
        }
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read as bracketed list

        // Slightly sub-optimal since it has intermediate resizing,
        // however don't expect this as input very often.

        list.clear();  // Clear addressing, leave storage intact (probably)

        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            // C++17
            // is >> list.emplace_back();

            // C++11
            list.emplace_back();
            is >> list.back();

            is.fatalCheck
            (
                "Istream >> std::vector<char> : "
                "reading entry"
            );

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }
    }
    else
    {
        list.clear();  // Clear old contents

        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const std::vector<T>& list)
{
    // Use UList for output
    UList<T> proxy(const_cast<T*>(list.data()), label(list.size()));
    os  << proxy;
    return os;
}


// ************************************************************************* //
