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

\*---------------------------------------------------------------------------*/

#include "PDRfitMeshParams.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRfitMeshParams::readBounds(const dictionary& dict)
{
    scalar value;

    if (dict.readIfPresent("xmin", value))
    {
        minBounds.x() = optionalData<scalar>(value);
    }
    if (dict.readIfPresent("ymin", value))
    {
        minBounds.y() = optionalData<scalar>(value);
    }
    if (dict.readIfPresent("zmin", value))
    {
        minBounds.z() = optionalData<scalar>(value);
    }

    if (dict.readIfPresent("xmax", value))
    {
        maxBounds.x() = optionalData<scalar>(value);
    }
    if (dict.readIfPresent("ymax", value))
    {
        maxBounds.y() = optionalData<scalar>(value);
    }
    if (dict.readIfPresent("zmax", value))
    {
        maxBounds.z() = optionalData<scalar>(value);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRfitMeshParams::PDRfitMeshParams(const dictionary& dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRfitMeshParams::read(const dictionary& dict)
{
    const dictionary* dictptr;

    // Geometric limits

    if ((dictptr = dict.findDict("bounds")) != nullptr)
    {
        readBounds(*dictptr);
    }

    if ((dictptr = dict.findDict("outer")) != nullptr)
    {
        const dictionary& d = *dictptr;

        d.readIfPresent("zmin", ground);
    }


    // Finding planes

    dict.readIfPresent("minFaceArea", minFaceArea);
    dict.readIfPresent("minAreaRatio", minAreaRatio);
    dict.readIfPresent("areaWidthFactor", areaWidthFactor);
    dict.readIfPresent("maxZoneToHeight", maxZoneToHeight);


    // Choosing cell width

    dict.readIfPresent("nCellsMin", nCellsMin);
    dict.readIfPresent("widthFactor", widthFactor);
    dict.readIfPresent("obsPerCell", obsPerCell);
    dict.readIfPresent("maxCellWidth", maxCellWidth);
    dict.readIfPresent("maxWidthEstimate", maxWidthEstimate);
    dict.readIfPresent("maxWidthRatio", maxWidthRatio);
    dict.readIfPresent("maxIterations", maxIterations);


    // Outer region

    dict.readIfPresent("nEdgeLayers", nEdgeLayers);

    dict.readIfPresent("outerRatio", outerRatio);

    // dict.readIfPresent("outerRadius", outerRadius);
};


// ************************************************************************* //
