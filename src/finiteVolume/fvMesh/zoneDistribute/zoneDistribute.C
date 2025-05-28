/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "zoneDistribute.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneDistribute, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneDistribute::zoneDistribute(const fvMesh& mesh)
:
    MeshObject_type(mesh),
    stencil_(zoneCPCStencil::New(mesh)),
    globalNumbering_(stencil_.globalNumbering()),
    send_(UPstream::nProcs()),
    pBufs_(UPstream::commsTypes::nonBlocking),
    cyclicBoundaryCells_(mesh.nCells(), false)
{
    // Don't clear storage on persistent buffer
    pBufs_.allowClearRecv(false);

    // Loop over boundary patches and store cells with a face on a cyclic patch
    bool hasCyclicPatches = false;
    forAll(mesh.boundaryMesh(), patchi)
    {
        const cyclicPolyPatch* cpp =
            isA<cyclicPolyPatch>(mesh.boundaryMesh()[patchi]);

        if (cpp)
        {
            cyclicBoundaryCells_.set(cpp->faceCells());
            hasCyclicPatches = true;
        }
    }

    // Populate cyclicCentres_
    if(hasCyclicPatches)
    {
        // Make a boolList from the bitSet
        boolList isCyclicCell(mesh.nCells(), false);

        forAll(cyclicBoundaryCells_, celli)
        {
            if (cyclicBoundaryCells_.test(celli))
            {
                isCyclicCell[celli] = true;
            }
        }

        // Use getFields to get map of cell centres across processor boundaries
        setUpCommforZone(isCyclicCell, true);

        cyclicCentres_.reset
        (
            new Map<vectorField>(getFields(isCyclicCell, mesh_.C()))
        );
    }
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::zoneDistribute& Foam::zoneDistribute::New(const fvMesh& mesh)
{
    auto* ptr = mesh.thisDb().getObjectPtr<zoneDistribute>("zoneDistribute");

    if (!ptr)
    {
        ptr = new zoneDistribute(mesh);
        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneDistribute::updateStencil(const boolList& zone)
{
    zoneCPCStencil::New(mesh_).updateStencil(zone);
}


void Foam::zoneDistribute::setUpCommforZone
(
    const boolList& zone,
    bool updateStencil
)
{
    zoneCPCStencil& stencil = zoneCPCStencil::New(mesh_);

    if (updateStencil)
    {
        stencil.updateStencil(zone);
    }

    if (UPstream::parRun())
    {
        List<labelHashSet> needed(UPstream::nProcs());

        // Bin according to originating (sending) processor
        for (const label celli : stencil.needsComm())
        {
            if (zone[celli])
            {
                for (const label gblIdx : stencil_[celli])
                {
                    const label proci = globalNumbering_.whichProcID(gblIdx);

                    if (proci != Pstream::myProcNo())
                    {
                        needed[proci].insert(gblIdx);
                    }
                }
            }
        }

        // Stream the send data into PstreamBuffers,
        // which we also use to track the current topology.

        pBufs_.clear();

        for (const int proci : pBufs_.allProcs())
        {
            const auto& indices = needed[proci];

            if (proci != UPstream::myProcNo() && !indices.empty())
            {
                // Serialize as List
                UOPstream toProc(proci, pBufs_);
                toProc << indices.sortedToc();
            }
        }

        pBufs_.finishedSends(sendConnections_, sendProcs_, recvProcs_);

        for (const int proci : pBufs_.allProcs())
        {
            send_[proci].clear();

            if (proci != UPstream::myProcNo() && pBufs_.recvDataCount(proci))
            {
                UIPstream fromProc(proci, pBufs_);
                fromProc >> send_[proci];
            }
        }
    }
}

Foam::List<Foam::label> Foam::zoneDistribute::getCyclicPatches
(
    const label celli,
    const label globalIdx,
    const vector globalIdxCellCentre
) const
{
    // Initialise cyclic patch label list
    List<label> patches(0);

    // If celli is not on a cyclic boundary, return the empty list
    if (!cyclicBoundaryCells_.test(celli))
    {
        return patches;
    }

    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    // Making list of cyclic patches to which celli belongs
    List<label> celliCyclicPatches;
    forAll(bMesh, patchi)
    {
        if (isA<cyclicPolyPatch>(bMesh[patchi]))
        {
            // Note: Probably not efficient due to use of found(celli) but
            // typically only used for very few cells (interface cells and their
            // point neighbours on cyclic boundaries).
            if (bMesh[patchi].faceCells().found(celli))
            {
                celliCyclicPatches.append(patchi);
            }
        }
    }

    // So celli belongs to at least one cyclic patch.
    // Let us figure out which.
    if (globalNumbering_.isLocal(globalIdx)) // celli and globalIdx on same proc
    {
        // Get all local point neighbor cells of celli, i.e. all point
        // neighbours that are not on the other side of a cyclic patch.
        List<label> localPointNeiCells(0);
        const labelList& cellPoints = mesh_.cellPoints()[celli];

        for (const label cellPoint : cellPoints)
        {
            const labelList& pointKCells = mesh_.pointCells()[cellPoint];

            for (const label pointKCell : pointKCells)
            {
                if (!localPointNeiCells.found(pointKCell))
                {
                    localPointNeiCells.append(pointKCell);
                }
            }
        }

        // Since globalIdx is a global cell index obtained from the point
        // neighbour list, stencil[celli], all cells in this that are not in
        // localPointNeiCells must be cyclic neighbour cells.
        const label localIdx = globalNumbering_.toLocal(globalIdx);
        if (!localPointNeiCells.found(localIdx))
        {
            for (const label patchi : celliCyclicPatches)
            {
                // Find the corresponding cyclic neighbor patch ID
                const cyclicPolyPatch& cpp =
                    static_cast<const cyclicPolyPatch&>(bMesh[patchi]);

                const label neiPatch = cpp.neighbPatchID();

                // Check if the cell globalIdx is on neiPatch.
                // If it is, append neiPatch to list of patches to return
                if (bMesh[neiPatch].faceCells().found(localIdx))
                {
                    patches.append(neiPatch);
                    // Here it may be possible to append patchi and do:
                    //
                    //    cpp.transformPosition()
                    //
                    // instead of
                    //
                    //    cpp.neighbPatch().transformPosition() in getPosition()
                }
            }
        }
    }
    else // celli and globalIdx on differet processors
    {
        List<label> cyclicID(3, -1);
        List<vector> separationVectors(3, vector(0,0,0));
        scalar distance = GREAT;

        forAll(celliCyclicPatches, cID)
        {
            cyclicID[cID] = celliCyclicPatches[cID];

            const label& patchI = celliCyclicPatches[cID];
            const cyclicPolyPatch& cpp =
                static_cast<const cyclicPolyPatch&>(bMesh[patchI]);

            if(cpp.transform() == coupledPolyPatch::transformType::ROTATIONAL)
            {
                FatalErrorInFunction
                    << "Rotational cyclic patches are not supported in parallel.\n"
                    << "Try to decompose the domain so that the rotational cyclic patch "
                    << "is not split in between processors."
                    << exit(FatalError);
            }
            cpp.neighbPatch().transformPosition(separationVectors[cID], 0);
        }

        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                for(int k = 0; k < 2; k++)
                {
                    vector separation =   i*separationVectors[0]
                                        + j*separationVectors[1]
                                        + k*separationVectors[2];

                    scalar testDistance = mag
                    (
                        (globalIdxCellCentre - separation)
                        -
                        mesh_.C()[celli]
                    );
                    if(debug) Info << "testDistance " << testDistance << endl;

                    if( testDistance < distance )
                    {
                        distance = testDistance;
                        patches = List<label>(0);
                        List<label> applyCyclic({i,j,k});

                        if(debug) Info << "distance " << distance << endl;
                        if(debug) Info << "separation " << separation << endl;
                        if(debug) Info << "applyCyclic " << applyCyclic << endl;

                        for(int n = 0; n < 3; n++)
                        {
                            if(cyclicID[n] != -1 && applyCyclic[n] == 1)
                            {
                                const cyclicPolyPatch& cpp =
                                    static_cast<const cyclicPolyPatch&>
                                    (
                                        bMesh[cyclicID[n]]
                                    );
                                if(debug)
                                {
                                    Info << "cpp.name() " << cpp.name() << endl;
                                }
                                patches.append(cpp.neighbPatchID());
                            }
                        }
                    }
                }
            }
        }
    }

    return patches;
}


Foam::vector Foam::zoneDistribute::getPosition
(
    const VolumeField<vector>& positions,
    const Map<vector>& valuesFromOtherProc,
    const label gblIdx,
    const List<label> cyclicPatchID
) const
{
    // Position vector, possibly from other processor, to be returned
    vector position(getValue(positions, valuesFromOtherProc, gblIdx));

    // Dealing with position transformation across cyclic patches.
    // If no transformation is required (most cases), cyclicPatchID is empty
    forAll(cyclicPatchID, i)
    {
        const label patchi = cyclicPatchID[i];

        const cyclicPolyPatch& cpp =
            static_cast<const cyclicPolyPatch&>
                (
                    positions.mesh().boundaryMesh()[patchi]
                );

        if (cpp.transform() != coupledPolyPatch::transformType::ROTATIONAL)
        {
            cpp.neighbPatch().transformPosition(position, 0);
        }
        else if (globalNumbering_.isLocal(gblIdx))
        {
            const label localIdx = globalNumbering_.toLocal(gblIdx);

            for (const label facei : mesh_.cells()[localIdx])
            {
                if (mesh_.boundaryMesh().whichPatch(facei) == cyclicPatchID[i])
                {
                    cpp.neighbPatch().transformPosition(position, facei);
                    continue;
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Rotational cyclic patches are not supported in parallel.\n"
                << "Try to decompose the domain so that the rotational cyclic"
                << "patch is not split in between processors."
                << exit(FatalError);
        }
    }

    return position;
}


// ************************************************************************* //
