/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
#include "DynamicList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcCellEdges() const
{
    // Loop through all faces and mark up cells with edges of the face.
    // Check for duplicates

    if (debug)
    {
        Pout<< "primitiveMesh::calcCellEdges() : "
            << "calculating cellEdges"
            << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // It is an error to attempt to recalculate cellEdges
    // if the pointer is already set
    if (cePtr_)
    {
        FatalErrorInFunction
            << "cellEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        // Set up temporary storage
        List<DynamicList<label>> ce(nCells());


        // Get reference to faceCells and faceEdges
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelListList& fe = faceEdges();

        // loop through the list again and add edges; checking for duplicates
        forAll(own, facei)
        {
            DynamicList<label>& curCellEdges = ce[own[facei]];

            const labelList& curEdges = fe[facei];

            for (const label edgei : curEdges)
            {
                // Add the edge
                if (!curCellEdges.contains(edgei))
                {
                    curCellEdges.push_back(edgei);
                }
            }
        }

        forAll(nei, facei)
        {
            DynamicList<label>& curCellEdges = ce[nei[facei]];

            const labelList& curEdges = fe[facei];

            for (const label edgei : curEdges)
            {
                // Add the edge
                if (!curCellEdges.contains(edgei))
                {
                    curCellEdges.push_back(edgei);
                }
            }
        }

        cePtr_ = std::make_unique<labelListList>(ce.size());
        auto& cellEdgeAddr = *cePtr_;

        // reset the size
        forAll(ce, celli)
        {
            cellEdgeAddr[celli].transfer(ce[celli]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::cellEdges() const
{
    if (!cePtr_)
    {
        calcCellEdges();
    }

    return *cePtr_;
}


// ************************************************************************* //
