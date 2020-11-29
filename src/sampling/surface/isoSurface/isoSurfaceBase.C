/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::isoSurfaceBase::algorithmType
>
Foam::isoSurfaceBase::algorithmNames
({
    { algorithmType::ALGO_CELL, "cell" },
    { algorithmType::ALGO_TOPO, "topo" },
    { algorithmType::ALGO_POINT, "point" },
});


const Foam::Enum
<
    Foam::isoSurfaceBase::filterType
>
Foam::isoSurfaceBase::filterNames
({
    { filterType::NONE, "none" },
    { filterType::CELL, "cell" },
    { filterType::DIAGCELL, "diagcell" },
    { filterType::PARTIAL, "partial" },
    { filterType::FULL, "full" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::isoSurfaceBase::algorithmType
Foam::isoSurfaceBase::getAlgorithmType
(
    const dictionary& dict,
    const isoSurfaceBase::algorithmType deflt
)
{
    word enumName;
    if (!dict.readIfPresent("isoAlgorithm", enumName, keyType::LITERAL))
    {
        return deflt;
    }

    if (!algorithmNames.found(enumName))
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << (algorithmNames) << nl
            << exit(FatalIOError);
    }

    return isoSurfaceBase::algorithmNames[enumName];
}


Foam::isoSurfaceBase::filterType
Foam::isoSurfaceBase::getFilterType
(
    const dictionary& dict,
    const isoSurfaceBase::filterType deflt
)
{
    word enumName;
    if (!dict.readIfPresent("regularise", enumName, keyType::LITERAL))
    {
        return deflt;
    }

    // Try as bool/switch
    const Switch sw = Switch::find(enumName);

    if (sw.good())
    {
        return (sw ? deflt : filterType::NONE);
    }

    // As enum
    if (!isoSurfaceBase::filterNames.found(enumName))
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << (filterNames) << nl
            << exit(FatalIOError);
    }

    return isoSurfaceBase::filterNames[enumName];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceBase::isoSurfaceBase
(
    const scalar iso,
    const boundBox& bounds
)
:
    meshedSurface(),
    iso_(iso),
    bounds_(bounds)
{}


// ************************************************************************* //
