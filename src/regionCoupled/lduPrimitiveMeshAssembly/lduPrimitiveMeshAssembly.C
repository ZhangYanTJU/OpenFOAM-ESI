/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "lduPrimitiveMeshAssembly.H"
#include "fvMesh.H"
#include "cyclicLduInterface.H"
#include "cyclicAssemblyFvPatch.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMeshAssembly, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::lduPrimitiveMeshAssembly::totalSize
(
    const UPtrList<fvMesh>& meshes
)
{
    label tot = 0;

    forAll(meshes, meshi)
    {
        tot += meshes[meshi].lduAddr().size();
    }

    return tot;
}


Foam::label Foam::lduPrimitiveMeshAssembly::findNbrMeshId
(
    const mappedPatchBase& pp,
    const UPtrList<fvMesh>& meshes
)
{
    forAll(meshes, meshi)
    {
        if (meshes[meshi].name() == pp.sampleMesh().name())
        {
            return meshi;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMeshAssembly::lduPrimitiveMeshAssembly
(
    const UPtrList<fvMesh>& meshes,
    const IOobject& io
)
:
    objectRegistry(io),
    lduPrimitiveMesh(totalSize(meshes))
{
    forAll(meshes, meshi)
    {
        if (meshes[meshi].comm() != comm())
        {
            WarningInFunction
                << "Communicator " << meshes[meshi].comm()
                << " at index " << meshi
                << " differs between meshes " << nl;
        }
    }

    const label nMeshes = meshes.size();
    patchMap_.setSize(nMeshes);
    faceBoundMap_.setSize(nMeshes);
    faceMap_.setSize(nMeshes);

    // Determine cellOffset and faceOffset
    cellOffsets_.setSize(1+nMeshes);
    cellOffsets_[0] = 0;
    for (label meshi=0; meshi < nMeshes; ++meshi)
    {
        cellOffsets_[meshi+1] =
            cellOffsets_[meshi] + meshes[meshi].lduAddr().size();
    }

    // Get newFaces and newPatches (not mapPolyPatches)
    label newFaces(0);
    label newPatches(0);
    for (label i=0; i < nMeshes; ++i)
    {
        patchMap_[i].setSize(meshes[i].boundaryMesh().size(), -1);

        forAll(meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (!isA<mappedWallPolyPatch>(pp))
            {
                patchMap_[i][patchI] = newPatches++;
            }
            else
            {
                const mappedPatchBase& mpp =
                    refCast<const mappedWallPolyPatch>(pp);

                const label meshNrbID = findNbrMeshId(mpp, meshes);

                const label nbrPatchID =
                    meshes[meshNrbID].boundaryMesh().findPatchID
                    (
                        mpp.samplePatch()
                    );

                const polyPatch& nbrpp =
                    meshes[meshNrbID].boundaryMesh()[nbrPatchID];

                if (pp.size() != nbrpp.size())
                {
                    FatalErrorInFunction
                        << "The number of faces on either side of the mapped"
                        << "patch " << pp.name() << " are not the same. "
                        << "This might be due to the decomposition used. "
                        << " Please use 'assembly' decomposition."
                        << exit(FatalError);
                }
                if (mpp.owner())
                {
                    newFaces += pp.size();
                }
            }
        }
    }

    // Add the internal faces for each mesh
    for (label i=0; i < nMeshes; ++i)
    {
        newFaces += meshes[i].lduAddr().upperAddr().size();
    }

    // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes[i].interfaces();
        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];
            if (globalPatchId != -1)
            {
                const labelUList& faceCells =
                    meshes[i].lduAddr().patchAddr(patchI);

                // Fill local patchAddr for standard patches
                if (!faceCells.empty())
                {
                    patchAddr_[globalPatchId].setSize(faceCells.size(), -1);

                    for (label celli = 0; celli < faceCells.size(); ++celli)
                    {
                        patchAddr_[globalPatchId][celli] =
                            cellOffsets_[i] + faceCells[celli];
                    }
                }
            }
        }
    }

    interfaces().setSize(newPatches);
    // Primitive interfaces
    primitiveInterfaces().setSize(newPatches);


    // The interfaces are conserved (cyclics, proc, etc)
    label interfaceID = 0;
    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes[i].interfaces();

        faceBoundMap_[i].setSize(interfacesLst.size());

        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];
            if (globalPatchId != -1)
            {
                // Set interfaces for cyclic (cyclicAssemblyFvPatch).
                // cyclic patches are cloned resetting nrbPatchID to the new
                // address plus local and nbr faceCells in global addressing
                // NOTE: the fvBoundary used remains the local. This means
                // that some member functions of the new cloned patch are
                // not fully functional, i.e neighbPatch() uses boundaryMesh
                // to access nbr patch.
                if (interfacesLst.set(patchI))
                {
                    if (isA<cyclicLduInterface>(interfacesLst[patchI]))
                    {
                        label nbrId = refCast
                            <const cyclicLduInterface>
                            (
                                interfacesLst[patchI]
                            ).neighbPatchID();

                        label globalNbr = patchMap()[i][nbrId];

                        primitiveInterfaces().set
                        (
                            interfaceID,
                            new cyclicAssemblyFvPatch
                            (
                               *(
                                    new cyclicPolyPatch
                                    (
                                        refCast<const cyclicPolyPatch>
                                        (
                                            meshes[i].boundaryMesh()[patchI]
                                        ),
                                        globalNbr,
                                        patchAddr_[globalPatchId]
                                    )
                                ),
                                meshes[i].boundary(),
                                patchAddr_[globalNbr]
                            )
                        );

                        interfaces().set
                        (
                            interfaceID,
                            &primitiveInterfaces()[interfaceID]
                        );
                    }
                    else
                    {
                        primitiveInterfaces().set
                        (
                            interfaceID,
                            nullptr
                        );

                        interfaces().set
                        (
                            interfaceID,
                            interfacesLst(patchI)
                        );
                    }
                }
                interfaceID++;
            }
        }
    }

    // Create new addressing
    lowerAddr().setSize(newFaces, -1);
    upperAddr().setSize(newFaces, -1);

    label startIndex = 0;

    for (label i=0; i < nMeshes; ++i)
    {
        faceMap_[i].setSize(meshes[i].lduAddr().lowerAddr().size(), -1);

        const label nFaces = meshes[i].lduAddr().upperAddr().size();

        // Add individual addresses
        SubList<label>(lowerAddr(), nFaces, startIndex) =
            meshes[i].lduAddr().lowerAddr();

        SubList<label>(upperAddr(), nFaces, startIndex) =
            meshes[i].lduAddr().upperAddr();

        // Offset cellsIDs to global cell addressing
        label localFacei = 0;

        for (label facei=startIndex; facei < startIndex + nFaces; ++facei)
        {
            lowerAddr()[facei] += cellOffsets_[i];
            upperAddr()[facei] += cellOffsets_[i];

            faceMap_[i][localFacei++] = facei;
        }

        startIndex += nFaces;
    }


    // Add lower/upper adressing for new internal faces corresponding
    // to old mapPolyPatch's
    label nFaces = startIndex;

    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (isA<mappedPatchBase>(pp))
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

                if (mpp.owner())
                {
                    labelList nbrFaceCells(mpp.samplePolyPatch().faceCells());

                    const label meshNrbId = findNbrMeshId(mpp, meshes);

                    if (meshNrbId == -1)
                    {
                         FatalErrorInFunction
                            << "Can not find Nbr Mesh for mapped patch"
                            << exit(FatalError);
                    }

                    forAll(pp.faceCells(), faceI)
                    {
                        const label cellI =
                            pp.faceCells()[faceI] + cellOffsets_[i];

                        const label nbrCellI =
                            nbrFaceCells[faceI] + cellOffsets_[meshNrbId];

                        lowerAddr()[nFaces] = min(cellI, nbrCellI);
                        upperAddr()[nFaces] = max(cellI, nbrCellI);

                        ++nFaces;
                    }
                }
            }
        }
    }

    if (newFaces != nFaces)
    {
       FatalErrorInFunction
            << "Incorrrect total number of faces in the assembled lduMatrix: "
            << newFaces << " != " << nFaces << nl
            << exit(FatalError);
    }

    // Fill faceBoundMap
    nFaces = startIndex;
    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (isA<mappedPatchBase>(pp) && !pp.empty())
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);
                if (mpp.owner())
                {
                    const label samplePatchi = mpp.samplePolyPatch().index();
                    const label meshNrbId = findNbrMeshId(mpp, meshes);
                    const polyPatch& ppNbr = mpp.samplePolyPatch();

                    if
                    (
                        faceBoundMap_[i][patchI].empty()
                     && faceBoundMap_[meshNrbId][samplePatchi].empty()
                    )
                    {
                        faceBoundMap_[i][patchI].setSize(pp.size(), -1);
                        faceBoundMap_[meshNrbId][samplePatchi].setSize
                        (
                            ppNbr.size(),
                            -1
                        );

                        forAll(pp.faceCells(), faceI)
                        {
                            faceBoundMap_[i][patchI][faceI] = nFaces;
                            faceBoundMap_[meshNrbId][samplePatchi][faceI] =
                                nFaces;
                            nFaces++;
                        }
                    }
                }
            }
        }
    }

    // Sort upper-tri order
    {
        labelList oldToNew
        (
            upperTriOrder
            (
                lduAddr().size(),
                lowerAddr(),
                upperAddr()
            )
        );

        inplaceReorder(oldToNew, lowerAddr());
        inplaceReorder(oldToNew, upperAddr());

        for (labelListList& bMap : faceBoundMap_)
        {
            for (labelList& faceMap : bMap)
            {
                for (label& facei : faceMap)
                {
                    facei = oldToNew[facei];
                }
            }
        }

        for (labelList& faceMap : faceMap_)
        {
            for (label& facei : faceMap)
            {
                facei = oldToNew[facei];
            }
        }
    }

    if (debug)
    {
        //DebugVar(faceBoundMap_);
        //DebugVar(lowerAddr());
        //DebugVar(upperAddr());
        //DebugVar(patchAddr_);
        //DebugVar(cellOffsets_);
        //DebugVar(faceMap_);
        checkUpperTriangular(lduAddr().size(), lowerAddr(), upperAddr());
    }
}


// ************************************************************************* //
