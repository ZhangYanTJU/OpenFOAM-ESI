/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class StringType, class UnaryPredicate>
StringType Foam::stringOps::quotemeta
(
    const StringType& str,
    const UnaryPredicate& meta,
    const char quote
)
{
    if (str.empty() || !quote)
    {
        return str;
    }

    StringType result;
    result.reserve(1.5*str.size());  // Moderately pessimistic

    bool escaped = false;
    for (const char c : str)
    {
        if (c == quote)
        {
            escaped = !escaped;  // toggle state
        }
        else if (escaped)
        {
            escaped = false;
        }
        else if (meta(c))
        {
            result += quote;
        }
        result += c;
    }

    result.shrink_to_fit();
    return result;
}


template<class StringType, class UnaryPredicate>
StringType Foam::stringOps::validate
(
    const std::string& str,
    const UnaryPredicate& accept,
    const bool invert
)
{
    StringType out;
    out.resize(str.length());

    std::string::size_type len = 0;

    for (std::string::size_type i = 0; i < str.length(); ++i)
    {
        const char c = str[i];
        if (accept(c) ? !invert : invert)
        {
            out[len++] += c;
        }
    }

    out.erase(len);
    return out;
}


// ************************************************************************* //
