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
#include "cyclicLduInterface.H"
#include "cyclicAssemblyFvPatch.H"
#include "mappedWallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "lduPrimitiveProcessorInterface.H"


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //


template<class Type>
void Foam::lduPrimitiveMeshAssembly::assemble()
{
    // Get newFaces and newPatches (not mapPolyPatches)
    label oldFaces(0);
    label newFaces(0);
    label newFacesProc(0);

    label nMeshes = meshes_.size();

    for (label i=0; i < nMeshes; ++i)
    {
        subFaceCompPatchMap_[i].setSize(meshes_[i].boundaryMesh().size());

        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];

            if (isA<mappedPatchBase>(pp))
            {
                const mappedPatchBase& mpp =
                    refCast<const mappedPatchBase>(pp);

                label meshNrbID = this->findNbrMeshId(mpp, meshes_);
                const label nbrPatchID =
                    meshes_[meshNrbID].boundaryMesh().findPatchID
                    (
                        mpp.samplePatch()
                    );

                const polyPatch& nbrpp =
                    meshes_[meshNrbID].boundaryMesh()[nbrPatchID];

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
            else if (isA<cyclicAMIPolyPatch>(pp))
            {
                const cyclicAMIPolyPatch& mpp =
                    refCast<const cyclicAMIPolyPatch>(pp);

                if (mpp.owner())
                {
                    const labelListList& addSourceFaces =
                        mpp.AMI().srcAddress();

                    List<DynamicList<label>> subFaceCompMap;

                    if (mpp.AMI().singlePatchProc() == -1)
                    {
                        FatalErrorInFunction
                            << "The number of faces on either side of the cyclicAMI"
                            << "patch " << pp.name() << " are not the same. "
                            << "This might be due to the decomposition used. "
                            << " Please use decomposition preserving the AMI on a "
                            << " single processor."
                            << exit(FatalError);

                        subFaceCompMap.setSize(mpp.AMI().tgtMap().constructSize());
                    }

                    // Add new faces as many weights for AMI
                    forAll (addSourceFaces, faceI)
                    {
                        const labelList& nbrFaceIs = addSourceFaces[faceI];

                        forAll (nbrFaceIs, j)
                        {
                            label nbrFaceI = nbrFaceIs[j];

                            if (nbrFaceI < mpp.neighbPatch().size())
                            {
                                // local faces
                                newFaces++;
                            }
                            else
                            {
                                if (mpp.AMI().singlePatchProc() == -1)
                                {
                                    // nbrFaceI is compactId paired with local faceI
                                    subFaceCompMap[nbrFaceI].append(faceI);
                                }
                                newFacesProc++;
                            }
                        }
                    }

                    subFaceCompPatchMap_[i][mpp.index()].setSize
                    (
                        subFaceCompMap.size()
                    );

                    labelListList& subFaceMap =
                        subFaceCompPatchMap_[i][mpp.index()];

                    forAll(subFaceMap, j)
                    {
                        subFaceMap[j] = std::move(subFaceCompMap[j]);
                    }
                }
            }
        }
    }


    // parallel. need to add procInterfaces between sub-faces and cells on
    // other proc's. Need to get local and remote faceCells to construct the
    // lduProcessorPrimitiveInterface
    reduce(newFacesProc, sumOp<label>());
    List<labelList> patchIntMap(meshes_.size());

    label interfaceI(0);
    label nbri(0);
    if (newFacesProc > 0)
    {
        List<labelList> nbrFaceCells(Pstream::nProcs());

        for (label i=0; i < nMeshes; ++i)
        {
            const polyBoundaryMesh& pbm = meshes_[i].boundaryMesh();

            patchIntMap[i].setSize(meshes_[i].boundaryMesh().size(), -1);
            myRmtTgtFaces_[i].setSize(meshes_[i].boundaryMesh().size());

            patchRemoteToLocal_[i].setSize(meshes_[i].boundaryMesh().size(), -1);

            forAll(pbm, patchi)
            {
                const polyPatch& pp = pbm[patchi];

                if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);

                    if (mpp.owner())
                    {
                        myRmtTgtFaces_[i][patchi].setSize(Pstream::nProcs());
                        myRmtTgtFaces_[i][mpp.neighbPatchID()].setSize(Pstream::nProcs());

                        List<DynamicList<label>> procOwner(Pstream::nProcs());
                        List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());

                        List<DynamicList<label>> dynTgtFaceIds(Pstream::nProcs());
                        List<DynamicList<label>> dynSrcFaceIds(Pstream::nProcs());

                        // Map per patchId per compactId for local faceIds
                        const labelListList& subFaceCompMap =
                            subFaceCompPatchMap_[i][mpp.index()];

                        const polyPatch& ppNbr = mpp.neighbPatch();
                        nbrFaceCells[Pstream::myProcNo()] = ppNbr.faceCells();

                        // Make available nbrFaceCells for all proc's
                        Pstream::gatherList(nbrFaceCells);
                        Pstream::scatterList(nbrFaceCells);

                        for (label proci = 0;proci<Pstream::nProcs(); proci++)
                        {
                            // Compact target map (key() = local in proci
                            // *iter = compactId)
                            const Map<label>& map = mpp.AMI().tgtcMap()[proci];
//DebugVar(map)
                            forAllConstIters(map, iter)
                            {
                                // Get nbrCell in proci
                                label cellI = nbrFaceCells[proci][iter.key()];

                                // Get local faces with stencil to other proc
                                const labelList& faces = subFaceCompMap[*iter];

                                forAll(faces, j)
                                {
                                    label localcellI =
                                        mpp.faceCells()[faces[j]];

                                    procOwner[proci].append(localcellI);
                                    dynProcNeighbour[proci].append(cellI);

                                    dynTgtFaceIds[proci].append(iter.key());
                                    dynSrcFaceIds[proci].append(faces[j]);
                                }
                            }
                        }
                        // Create proc interface per cyclicAMI if needed
                        createRemoteInterface
                        (
                            procOwner,
                            dynProcNeighbour,
                            patchRemoteToLocal_[i],
                            mpp.index(),
                            interfaceI
                        );

                        labelListList dynTgtFace(dynTgtFaceIds.size());
                        forAll(dynTgtFace, j)
                        {
                            dynTgtFace[j] = std::move(dynTgtFaceIds[j]);
                        }

                        labelListList mySrcFace(dynSrcFaceIds.size());
                        forAll(mySrcFace, j)
                        {
                            mySrcFace[j] = std::move(dynSrcFaceIds[j]);
                        }

                        forAll(mySrcFace, j)
                        {
                            myRmtTgtFaces_[i][mpp.index()][j].setSize
                            (
                                mySrcFace[j].size()
                            );
                            myRmtTgtFaces_[i][mpp.index()][j] = std::move(mySrcFace[j]);
                        }

                        // myRmtTgtFaces_ holds slave faces on myProc
                        Pstream::exchange<labelList, label>
                        (
                            dynTgtFace,
                            myRmtTgtFaces_[i][mpp.neighbPatchID()]
                        );
                    }
                }
            }
        }
        remoteStencilInterfaces_.shrink();
        nbri = remoteStencilInterfaces_.size();
    }

    label newPatches(0);

    // patchMap: maps from original to asembled (-1 for removed)
    // patchLocalToGlobalMap: map from original to asembled + extra Ids
    // for flux treatment on original Ids.
    // patchRemoteToLocal: maps from remoteI to removed patchId. Needed
    // for setting internal and boundary coeffs.
    for (label i=0; i < nMeshes; ++i)
    {
        patchMap_[i].setSize(meshes_[i].boundaryMesh().size() + nbri, -1);
        patchLocalToGlobalMap_[i].setSize(patchMap_[i].size(), -1);

        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];

            if (!isA<Type>(pp))
            {
                patchMap_[i][patchI] = newPatches;
                patchLocalToGlobalMap_[i][patchI] = newPatches;
                newPatches++;
            }
        }

        label firstRemote = meshes_[i].boundaryMesh().size();
        forAll(remoteStencilInterfaces_, remoteI)
        {
            if (remoteStencilInterfaces_.set(remoteI))
            {
                patchMap_[i][firstRemote] = newPatches;
                patchLocalToGlobalMap_[i][firstRemote] = newPatches;
                newPatches++;
                firstRemote++;
            }
        }
    }

    label virtualPatches = newPatches;
    //DebugVar(patchMap_)
    //DebugVar(patchRemoteToLocal_)

    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            if (patchLocalToGlobalMap_[i][patchI] == -1)
            {
                patchLocalToGlobalMap_[i][patchI] = virtualPatches++;
            }
        }
    }

    //DebugVar(patchLocalToGlobalMap_)

    // Add the internal faces for each mesh
    for (label i=0; i < nMeshes; ++i)
    {
        newFaces += meshes_[i].lduAddr().upperAddr().size();
        oldFaces += meshes_[i].lduAddr().upperAddr().size();
    }

    if (debug)
    {
        Pout<< "old nFaces : " << oldFaces
            << "new nFaces (internal): " << newFaces
            << "new nFaces (remote) : " << newFacesProc << endl;

        Pout<< "new patches : " << newPatches << endl;
    }
    // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes_[i].interfaces();

        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];

            if (globalPatchId != -1)
            {
                const labelUList& faceCells =
                    meshes_[i].lduAddr().patchAddr(patchI);

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

        // Set patchAddr for remote interfaces
        label firstRemote = meshes_[i].boundaryMesh().size();

        for (label j=0; j<nbri; j++)
        {
            label patchI = patchMap_[i][firstRemote];

            const lduPrimitiveProcessorInterface& pp =
                remoteStencilInterfaces_[j];

            const labelUList& faceCells = pp.faceCells();

            patchAddr_[patchI].setSize(faceCells.size(), -1);

            for (label celli = 0; celli < faceCells.size(); ++celli)
            {
                patchAddr_[patchI][celli] = faceCells[celli];
            }
            firstRemote++;

            //Pout<< "cells remote : " << patchI << " "
            //    << patchAddr_[patchI] << endl;
        }
    }

    // Interfaces
    interfaces().setSize(newPatches);
    // Primitive interfaces
    primitiveInterfaces().setSize(newPatches);

    // The interfaces are conserved (cyclics, proc, etc)
    label interfaceID = 0;
    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes_[i].interfaces();

        faceBoundMap_[i].setSize(interfacesLst.size());
        cellBoundMap_[i].setSize(interfacesLst.size());
        magSfFaceBoundMap_[i].setSize(interfacesLst.size());

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
                    if
                    (
                        isA<cyclicLduInterface>(interfacesLst[patchI])
                     && resetCyclics_
                    )
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
                                            meshes_[i].boundaryMesh()[patchI]
                                        ),
                                        globalNbr,
                                        patchAddr_[globalPatchId]
                                    )
                                ),
                                meshes_[i].boundary(),
                                patchAddr_[globalNbr]
                            )
                        );

                        interfaces().set
                        (
                            interfaceID,
                            &primitiveInterfaces()[interfaceID]
                        );
                    }
                    else // All the other patches are conserved
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

        forAll(remoteStencilInterfaces_, remoteI)
        {
            if (remoteStencilInterfaces_.set(remoteI))
            {
                const lduPrimitiveProcessorInterface& pp =
                    remoteStencilInterfaces_[remoteI];

                interfaces().set(interfaceID, &pp);
                primitiveInterfaces().set
                (
                    interfaceID,
                    nullptr
                );
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
        faceMap_[i].setSize(meshes_[i].lduAddr().lowerAddr().size(), -1);

        const label nFaces = meshes_[i].lduAddr().upperAddr().size();

        // Add individual addresses
        SubList<label>(lowerAddr(), nFaces, startIndex) =
            meshes_[i].lduAddr().lowerAddr();

        SubList<label>(upperAddr(), nFaces, startIndex) =
            meshes_[i].lduAddr().upperAddr();

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


    // Add new lower/upper adressing for new internal faces corresponding
    // to patch faces that has a correspondent on the slave patch
    // (i.e map, AMI,etc)
    // Don't include faces that are in different proc
    label nFaces = startIndex;

    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];

            if (isA<Type>(pp))
            {
                label meshNrbId = 0;
                labelList nbrFaceCells;

                if (isA<mappedPatchBase>(pp))
                {
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>(pp);

                    if (mpp.owner())
                    {
                        cellBoundMap_[i][patchI].setSize(pp.size(), -1);

                        nbrFaceCells =
                        (
                            mpp.samplePolyPatch().faceCells()
                        );

                        const label nbrPatchId =
                            mpp.samplePolyPatch().index();

                        cellBoundMap_[i][nbrPatchId].setSize(pp.size(), -1);

                        meshNrbId = this->findNbrMeshId(mpp, meshes_);

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

                            cellBoundMap_[i][patchI][faceI] = nbrCellI;
                            cellBoundMap_[i][nbrPatchId][faceI] = cellI;

                            ++nFaces;
                        }
                    }
                }
                else if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);

                    if (mpp.owner() && mpp.size() > 0)
                    {
                        const labelListList& sourceFaces =
                            mpp.AMI().srcAddress();

                        nbrFaceCells =
                        (
                            mpp.neighbPatch().faceCells()
                        );

                        const label nbrPatchId = mpp.neighbPatchID();

                        label newFaces(0);
                        label pFaces(0);
                        forAll (sourceFaces, faceI)
                        {
                            const labelList& faceIds = sourceFaces[faceI];
                            forAll (faceIds, j)
                            {
                                if (faceIds[j] < mpp.neighbPatch().size())
                                {
                                    newFaces++;
                                }
                                else
                                {
                                    pFaces++;
                                }
                            }
                        }

                        cellBoundMap_[i][patchI].setSize(newFaces + pFaces, -1);
                        cellBoundMap_[i][nbrPatchId].setSize(newFaces + pFaces, -1);

                        magSfFaceBoundMap_[i][patchI].setSize(newFaces + pFaces, -1);
                        magSfFaceBoundMap_[i][nbrPatchId].setSize
                        (
                            newFaces + pFaces, -1
                        );

                        // Compact target map
                        // key() = local face in proci
                        // *iter = compactId

                        label subFaceI(0);
                        forAll(pp.faceCells(), faceI)
                        {
                            const label cellI =
                                pp.faceCells()[faceI] + cellOffsets_[i];

                            const labelList& facesIds = sourceFaces[faceI];

                            forAll(facesIds, j)
                            {
                                label nbrFaceId = facesIds[j];

                                if (nbrFaceId < mpp.neighbPatch().size())
                                {
                                    // local faces
                                    const label nbrCellI  = nbrFaceCells[nbrFaceId]
                                        + cellOffsets_[meshNrbId];

                                    lowerAddr()[nFaces] = min(cellI, nbrCellI);
                                    upperAddr()[nFaces] = max(cellI, nbrCellI);

                                    cellBoundMap_[i][patchI][subFaceI] = nbrCellI;
                                    cellBoundMap_[i][nbrPatchId][subFaceI] = cellI;

                                    magSfFaceBoundMap_[i][patchI][subFaceI] =
                                        faceI;
                                    magSfFaceBoundMap_[i][nbrPatchId][subFaceI] =
                                        nbrFaceId;
                                    //faceI;

                                    ++subFaceI;
                                    ++nFaces;
                                }
                            }
                        }
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
        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];
            if (isA<Type>(pp))
            {
                label meshNrbId = 0;
                label samplePatchi = -1;

                if (isA<mappedPatchBase>(pp))
                {
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>(pp);
                    if (mpp.owner())
                    {
                        samplePatchi = mpp.samplePolyPatch().index();
                        meshNrbId = this->findNbrMeshId(mpp, meshes_);

                        if
                        (
                            faceBoundMap_[i][patchI].empty()
                         && faceBoundMap_[meshNrbId][samplePatchi].empty()
                        )
                        {
                            faceBoundMap_[i][patchI].setSize(pp.size(), -1);
                            faceBoundMap_[meshNrbId][samplePatchi].setSize
                                (pp.size(),-1);
                        }

                        forAll(pp.faceCells(), faceI)
                        {
                            faceBoundMap_[i][patchI][faceI] = nFaces;
                            faceBoundMap_[meshNrbId][samplePatchi][faceI] =
                                nFaces;
                            nFaces++;
                        }
                    }
                }
                else if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);
                    if (mpp.owner())
                    {
                        samplePatchi = mpp.neighbPatchID();

                        if (faceBoundMap_[i][patchI].empty())
                        {
                            faceBoundMap_[i][patchI].setSize
                            (
                                cellBoundMap_[i][patchI].size(),
                                -1
                            );

                            faceBoundMap_[i][samplePatchi].setSize
                            (
                                cellBoundMap_[i][samplePatchi].size(),
                                -1
                            );

                            const labelListList& sourceFaces =
                                mpp.AMI().srcAddress();

                            label subFaceI(0);
                            forAll(pp.faceCells(), faceI)
                            {
                                const labelList& facesIds = sourceFaces[faceI];

                                // nFaces are the new sub-facesIds in the new global index
                                // subFaceI are the local sub-faceIds
                                forAll(facesIds, j)
                                {
                                    if (facesIds[j] <  mpp.neighbPatch().size())
                                    {
                                        faceBoundMap_[i][patchI][subFaceI] =
                                            nFaces;
                                        faceBoundMap_[i][samplePatchi][subFaceI] =
                                            nFaces;
                                        subFaceI++;
                                        nFaces++;
                                    }
                                }

                            }
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

        for (labelList& faceMap : faceMap_)
        {
            for (label& facei : faceMap)
            {
                facei = oldToNew[facei];
            }
        }

        for (labelListList& bMap : faceBoundMap_)
        {
            for (labelList& faceMap : bMap)
            {
                for (label& facei : faceMap)
                {
                    if (facei != -1)
                    {
                        facei = oldToNew[facei];
                    }
                }
            }
        }
    }

    if (debug)
    {
        DebugVar(faceBoundMap_);
        DebugVar(lowerAddr());
        DebugVar(upperAddr());
        DebugVar(patchAddr_);
        DebugVar(cellOffsets_);
        DebugVar(faceMap_);
        checkUpperTriangular(lduAddr().size(), lowerAddr(), upperAddr());
    }
}
// ************************************************************************* //
