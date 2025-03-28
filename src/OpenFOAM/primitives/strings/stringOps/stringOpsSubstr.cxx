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

\*---------------------------------------------------------------------------*/

#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::SubStrings Foam::stringOps::split
(
    const std::string& str,
    const char delim,
    std::string::size_type pos,
    const bool keepEmpty
)
{
    Foam::SubStrings list;

    if
    (
        !delim
     || (pos == std::string::npos || pos >= str.size())
    )
    {
        return list;
    }

    list.reserve(20);

    std::string::size_type end;
    while ((end = str.find(delim, pos)) != std::string::npos)
    {
        if (keepEmpty || (pos < end))
        {
            list.append(str.cbegin() + pos, str.cbegin() + end);
        }
        pos = end + 1;
    }

    // Trailing element
    if (keepEmpty ? (pos <= str.size()) : (pos < str.size()))
    {
        list.append(str.cbegin() + pos, str.cend());
    }

    return list;
}


Foam::SubStrings Foam::stringOps::split
(
    const std::string& str,
    const std::string& delim,
    std::string::size_type pos,
    const bool keepEmpty
)
{
    Foam::SubStrings list;

    if
    (
        delim.empty()
     || (pos == std::string::npos || pos >= str.size())
    )
    {
        return list;
    }

    list.reserve(20);

    std::string::size_type end;
    while ((end = str.find(delim, pos)) != std::string::npos)
    {
        if (keepEmpty || (pos < end))
        {
            list.append(str.cbegin() + pos, str.cbegin() + end);
        }
        pos = end + delim.size();
    }

    // Trailing element
    if (keepEmpty ? (pos <= str.size()) : (pos < str.size()))
    {
        list.append(str.cbegin() + pos, str.cend());
    }

    return list;
}


Foam::SubStrings Foam::stringOps::splitAny
(
    const std::string& str,
    const std::string& delim,
    std::string::size_type pos
)
{
    Foam::SubStrings list;

    if
    (
        delim.empty()
     || (pos == std::string::npos || pos >= str.size())
    )
    {
        return list;
    }

    list.reserve(20);

    while ((pos = str.find_first_not_of(delim, pos)) != std::string::npos)
    {
        const auto end = str.find_first_of(delim, pos);

        if (end == std::string::npos)
        {
            // Trailing element
            list.append(str.cbegin() + pos, str.cend());
            break;
        }

        // Intermediate element
        list.append(str.cbegin() + pos, str.cbegin() + end);

        pos = end + 1;
    }

    return list;
}


Foam::SubStrings Foam::stringOps::splitFixed
(
    const std::string& str,
    const std::string::size_type width,
    std::string::size_type pos
)
{
    Foam::SubStrings list;

    if
    (
        !width
     || (pos == std::string::npos || pos >= str.size())
    )
    {
        return list;
    }

    list.reserve(1 + ((str.size() - pos) / width));

    const auto len = str.size();

    while (pos < len)
    {
        const auto end = (pos + width);

        if (end >= len)
        {
            // Trailing element
            list.append(str.cbegin() + pos, str.cend());
            break;
        }

        // Intermediate element
        list.append(str.cbegin() + pos, str.cbegin() + end);

        pos += width;
    }

    return list;
}


Foam::SubStrings Foam::stringOps::splitSpace
(
    const std::string& str,
    std::string::size_type pos
)
{
    return splitAny(str, "\t\n\v\f\r ", pos);
}


// ************************************************************************* //
