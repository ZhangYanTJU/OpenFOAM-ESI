/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
#include "cell.H"
#include "bitSet.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcPointCells() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcPointCells() : "
            << "calculating pointCells"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate pointCells
    // if the pointer is already set
    if (pcPtr_)
    {
        FatalErrorInFunction
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else
    {
        // Invert cellPoints
        (void)cellPoints();
        pcPtr_ = new labelListList(nPoints());
        invertManyToMany(nPoints(), cellPoints(), *pcPtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::pointCells() const
{
    if (!pcPtr_)
    {
        calcPointCells();
    }

    return *pcPtr_;
}


const Foam::labelList& Foam::primitiveMesh::pointCells
(
    const label pointi,
    DynamicList<label>& storage
) const
{
    if (hasPointCells())
    {
        return pointCells()[pointi];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelList& pFaces = pointFaces()[pointi];

        storage.clear();

        for (const label facei : pFaces)
        {
            // Owner cell
            storage.push_back(own[facei]);

            // Neighbour cell
            if (facei < nInternalFaces())
            {
                storage.push_back(nei[facei]);
            }
        }

        // Filter duplicates
        if (storage.size() > 1)
        {
            std::sort(storage.begin(), storage.end());

            auto last = std::unique(storage.begin(), storage.end());

            const label newLen = label(last - storage.begin());

            storage.resize(newLen);
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::pointCells(const label pointi) const
{
    return pointCells(pointi, labels_);
}


// ************************************************************************* //
