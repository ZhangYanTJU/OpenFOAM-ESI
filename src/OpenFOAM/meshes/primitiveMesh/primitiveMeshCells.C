/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCells
(
    cellList& cellFaceAddr,
    const labelUList& own,
    const labelUList& nei,
    const label inNCells
)
{
    label nCells = inNCells;

    if (nCells == -1)
    {
        nCells = -1;

        forAll(own, facei)
        {
            nCells = max(nCells, own[facei]);
        }
        nCells++;
    }

    // 1. Count number of faces per cell

    labelList ncf(nCells, Zero);

    forAll(own, facei)
    {
        ncf[own[facei]]++;
    }

    forAll(nei, facei)
    {
        if (nei[facei] >= 0)
        {
            ncf[nei[facei]]++;
        }
    }

    // Create the storage
    cellFaceAddr.setSize(ncf.size());


    // 2. Size and fill cellFaceAddr

    forAll(cellFaceAddr, celli)
    {
        cellFaceAddr[celli].setSize(ncf[celli]);
    }
    ncf = 0;

    forAll(own, facei)
    {
        label celli = own[facei];

        cellFaceAddr[celli][ncf[celli]++] = facei;
    }

    forAll(nei, facei)
    {
        label celli = nei[facei];

        if (celli >= 0)
        {
            cellFaceAddr[celli][ncf[celli]++] = facei;
        }
    }
}


void Foam::primitiveMesh::calcCells() const
{
    // Loop through faceCells and mark up neighbours

    if (debug)
    {
        Pout<< "primitiveMesh::calcCells() : calculating cells"
            << endl;
    }

    // It is an error to attempt to recalculate cells
    // if the pointer is already set
    if (cfPtr_)
    {
        FatalErrorInFunction
            << "cells already calculated"
            << abort(FatalError);
    }
    else
    {
        // Create the storage
        cfPtr_ = std::make_unique<cellList>(nCells());
        auto& cellFaceAddr = *cfPtr_;

        calcCells
        (
            cellFaceAddr,
            faceOwner(),
            faceNeighbour(),
            nCells()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cellList& Foam::primitiveMesh::cells() const
{
    if (!cfPtr_)
    {
        calcCells();
    }

    return *cfPtr_;
}


// ************************************************************************* //
