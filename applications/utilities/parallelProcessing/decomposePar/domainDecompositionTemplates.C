/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BinaryOp>
void Foam::domainDecomposition::processInterCyclics
(
    const polyBoundaryMesh& patches,
    List<DynamicList<DynamicList<label>>>& interPatchFaces,
    List<Map<label>>& procNbrToInterPatch,
    List<labelListList>& subPatchIDs,
    List<labelListList>& subPatchStarts,
    bool owner,
    BinaryOp bop
) const
{
    // Processor boundaries from split cyclics
    forAll(patches, patchi)
    {
        const auto& pp = patches[patchi];
        const auto* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner() == owner)
        {
            // cyclic: check opposite side on this processor
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            // cyclic: check opposite side on this processor
            const labelUList& patchFaceCells = pp.faceCells();
            const labelUList& nbrPatchFaceCells = nbrPatch.faceCells();

            // Store old sizes. Used to detect which inter-proc patches
            // have been added to.
            labelListList oldInterfaceSizes(nProcs_);
            forAll(oldInterfaceSizes, proci)
            {
                labelList& curOldSizes = oldInterfaceSizes[proci];

                curOldSizes.setSize(interPatchFaces[proci].size());
                forAll(curOldSizes, interI)
                {
                    curOldSizes[interI] =
                        interPatchFaces[proci][interI].size();
                }
            }

            // Add faces with different owner and neighbour processors
            forAll(patchFaceCells, facei)
            {
                const label ownerProc = cellToProc_[patchFaceCells[facei]];
                const label nbrProc = cellToProc_[nbrPatchFaceCells[facei]];
                if (bop(ownerProc, nbrProc))
                {
                    // inter - processor patch face found.
                    addInterProcFace
                    (
                        pp.start()+facei,
                        ownerProc,
                        nbrProc,
                        procNbrToInterPatch,
                        interPatchFaces
                    );
                }
            }

            // 1. Check if any faces added to existing interfaces
            forAll(oldInterfaceSizes, proci)
            {
                const labelList& curOldSizes = oldInterfaceSizes[proci];

                forAll(curOldSizes, interI)
                {
                    label oldSz = curOldSizes[interI];
                    if (interPatchFaces[proci][interI].size() > oldSz)
                    {
                        // Added faces to this interface. Add an entry
                        subPatchIDs[proci][interI].append(patchi);
                        subPatchStarts[proci][interI].append(oldSz);
                    }
                }
            }

            // 2. Any new interfaces
            forAll(subPatchIDs, proci)
            {
                label nIntfcs = interPatchFaces[proci].size();
                subPatchIDs[proci].setSize(nIntfcs, labelList(1, patchi));
                subPatchStarts[proci].setSize(nIntfcs, labelList(1, Zero));
            }
        }
    }
}


// ************************************************************************* //
