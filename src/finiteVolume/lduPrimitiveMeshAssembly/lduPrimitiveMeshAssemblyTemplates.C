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

    labelListListList subFaceCompPatchMap(meshes_[0].boundaryMesh().size());

    for (label i=0; i < nMeshes; ++i)
    {
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

                    List<DynamicList<label>>
                        subFaceCompMap(mpp.AMI().tgtMap().constructSize());

                    // Add new faces as many weights for AMI
                    forAll (addSourceFaces, faceI)
                    {
                        const labelList& nbrFaceIs = addSourceFaces[faceI];

                        forAll (nbrFaceIs, i)
                        {
                            label nbrFaceI = nbrFaceIs[i];

                            if (nbrFaceI < mpp.neighbPatch().size())
                            {
                                // local faces
                                newFaces++;
                            }
                            else
                            {   //nbrFaceI is compactId map with local faceI
                                subFaceCompMap[nbrFaceI].append(faceI);
                                newFacesProc++;
                            }
                        }

                        Pout << "faceI " << faceI << " "
                            << mpp.AMI().srcAddress()[faceI] << endl;
                    }

                    subFaceCompPatchMap[mpp.index()].setSize(subFaceCompMap.size());
                    labelListList& subFaceMap = subFaceCompPatchMap[mpp.index()];
                    forAll(subFaceMap, i)
                    {
                        subFaceMap[i] = std::move(subFaceCompMap[i]);
                    }
                }
            }
        }
    }

    // Per processor to owner (local)/neighbour (remote)
    List<DynamicList<label>> procOwner(Pstream::nProcs());
    List<DynamicList<label>> dynProcNeighbour(Pstream::nProcs());

    // parallel. need to add procInterfaces between sub-faces and cells on other
    // proc's. Need to get faceCells to construct the lduProcessorPrimitiveInterface
    if (newFacesProc > 0)
    {
        List<labelList> nbrFaceCells(Pstream::nProcs());

        for (label i=0; i < nMeshes; ++i)
        {
            const polyBoundaryMesh& pbm = meshes_[i].boundaryMesh();
            forAll(pbm, patchi)
            {
                const polyPatch& pp = pbm[patchi];

                if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);

                    if (mpp.owner())
                    {
                        const labelListList& subFaceCompMap =
                            subFaceCompPatchMap[mpp.index()];

                        const polyPatch& ppNbr = mpp.neighbPatch();
                        nbrFaceCells[Pstream::myProcNo()] = ppNbr.faceCells();
DebugVar(nbrFaceCells)
                        Pstream::gatherList(nbrFaceCells);
                        Pstream::scatterList(nbrFaceCells);
DebugVar(nbrFaceCells)
                        for (label proci = 0;proci<Pstream::nProcs(); proci++)
                        {
                            //if (proci != Pstream::myProcNo())
                            {
                                const Map<label>& map =
                                    mpp.AMI().tgtcMap()[proci];
DebugVar(map)
                                forAllConstIters(map, iter)
                                {
DebugVar(iter.key())
                                    label cellI = nbrFaceCells[proci][iter.key()];
DebugVar(*iter)
                                    const labelList& faces =
                                        subFaceCompMap[*iter];
DebugVar(faces)
                                    forAll(faces, faceI)
                                    {
                                        label localcellI =
                                            mpp.faceCells()[faces[faceI]];
                                        procOwner[proci].append
                                        (
                                            localcellI
                                        );
                                        dynProcNeighbour[proci].append(cellI);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    labelListList procNeighbour(dynProcNeighbour.size());
    forAll(procNeighbour, i)
    {
        procNeighbour[i] = std::move(dynProcNeighbour[i]);
    }
    DebugVar(procNeighbour)

    labelListList mySendCells;
    Pstream::exchange<labelList, label>(procNeighbour, mySendCells);

    DebugVar(mySendCells)

    DebugVar(procOwner)

    label nbri = 0;
    forAll(procOwner, proci)
    {
        if (procOwner[proci].size())
        {
            nbri++;
        }
        if (mySendCells[proci].size())
        {
            nbri++;
        }
    }
    remoteStencilInterfaces_.setSize(nbri);

    nbri = 0;

    labelList procToInterface(Pstream::nProcs(), -1);
    forAll(procOwner, proci)
    {
        if (proci < Pstream::myProcNo() && procOwner[proci].size())
        {
            if (1)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci]
                    << " from " << proci
                    << " tag " << Pstream::msgType()+2 << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
        else if (proci > Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (1)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci]
                    << " to " << proci
                    << " tag " << Pstream::msgType()+2 << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+2
                )
            );
        }
    }
    forAll(procOwner, proci)
    {
        if (proci > Pstream::myProcNo() && procOwner[proci].size())
        {
            if (1)
            {
                Pout<< "Adding interface " << nbri
                    << " to receive my " << procOwner[proci]
                    << " from " << proci
                    << " tag " << Pstream::msgType()+3 << endl;
            }
            procToInterface[proci] = nbri;
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    procOwner[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
        else if (proci < Pstream::myProcNo() && mySendCells[proci].size())
        {
            if (1)
            {
                Pout<< "Adding interface " << nbri
                    << " to send my " << mySendCells[proci]
                    << " to " << proci
                    << " tag " << Pstream::msgType()+3 << endl;
            }
            remoteStencilInterfaces_.set
            (
                nbri++,
                new lduPrimitiveProcessorInterface
                (
                    mySendCells[proci],
                    Pstream::myProcNo(),
                    proci,
                    tensorField(0),
                    Pstream::msgType()+3
                )
            );
        }
    }

    label newPatches(0);

    // Size patchMap amd patchLocalToGlobalMap
    // Need to find out which type of patch generated the remote
    // to be called to fill boundary/internal coeffs
    for (label i=0; i < nMeshes; ++i)
    {
        patchMap_[i].setSize(meshes_[i].boundaryMesh().size() + nbri, -1);
        patchLocalToGlobalMap_[i].setSize(patchMap_[i].size(), -1);
        patchRemoteToLocal_[i].setSize(remoteStencilInterfaces_.size(), -1);

        label nRemoteI = 0;
        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];

            if (!isA<Type>(pp))
            {
                patchMap_[i][patchI] = newPatches;
                patchLocalToGlobalMap_[i][patchI] = newPatches;
                newPatches++;
            }
            else if (isA<cyclicAMIPolyPatch>(pp))
            {
                const cyclicAMIPolyPatch& mpp =
                    refCast<const cyclicAMIPolyPatch>(pp);

                if (mpp.owner())
                {
                    patchRemoteToLocal_[i][nRemoteI++] = mpp.index();
                }
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
    DebugVar(patchMap_)
    DebugVar(patchRemoteToLocal_)

    // patchLocalToGlobalMap local -> removed patches Id's
    // This is used for the flux on patches which were removed
    // but used in the .flux() of the original matrix
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

    DebugVar(patchLocalToGlobalMap_)

    // Add the internal faces for each mesh
    for (label i=0; i < nMeshes; ++i)
    {
        newFaces += meshes_[i].lduAddr().upperAddr().size();
        oldFaces += meshes_[i].lduAddr().upperAddr().size();
    }

    Pout<< "old nFaces : " << oldFaces
        << "new nFaces (internal): " << newFaces
        << "new nFaces (remote) : " << newFacesProc << endl;

    // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    Pout<< "new patches : " << newPatches << endl;

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

        label firstRemote = meshes_[i].boundaryMesh().size();
        //forAll(remoteStencilInterfaces_, remoteI)
        for (label j=0; j<nbri; j++)
        {
            label patchI = patchMap_[i][firstRemote];

            const lduPrimitiveProcessorInterface& pp = remoteStencilInterfaces_[j];

            const labelUList& faceCells = pp.faceCells();

            patchAddr_[patchI].setSize(faceCells.size(), -1);

            for (label celli = 0; celli < faceCells.size(); ++celli)
            {
                patchAddr_[patchI][celli] = faceCells[celli];
            }
            firstRemote++;
            DebugVar(patchAddr_[patchI])
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
    // to patch faces that has a correspondent on the slave patch (i.e map, AMI,etc)
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

                    if (mpp.owner())
                    {
                        const labelListList& sourceFaces =
                            mpp.AMI().srcAddress();

                        nbrFaceCells =
                        (
                            mpp.neighbPatch().faceCells()
                        );

                        const label nbrPatchId = mpp.neighbPatchID();

                        label newFaces(0);
                        label procFaces(0);
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
                                    procFaces++;
                                }
                            }
                                //newFaces += sourceFaces[faceI].size();
                        }

DebugVar(newFaces)
DebugVar(procFaces)
                        cellBoundMap_[i][patchI].setSize(newFaces + procFaces, -1);
                        cellBoundMap_[i][nbrPatchId].setSize(newFaces + procFaces, -1);
                        magSfFaceBoundMap_[i][patchI].setSize(newFaces + procFaces, -1);
                        magSfFaceBoundMap_[i][nbrPatchId].setSize(newFaces + procFaces, -1);

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
                                else
                                {
                                    Pout << "face " << nbrFaceId << " other proc" << endl;
                                }
                            }
                        }
        //Pout << "cellBoundMap patch " << cellBoundMap_[i][patchI] << endl;
        //Pout << "cellBoundMap nbrpatch " << cellBoundMap_[i][nbrPatchId] << endl;
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
//Pout << "faceBoundMap_ " <<  faceBoundMap_[i][patchI] << endl;
                    }
//                     else
//                     {
//                         if (faceBoundMap_[i][patchI].empty())
//                         {
//                             faceBoundMap_[i][patchI].setSize
//                             (
//                                 cellBoundMap_[i][patchI].size(),
//                                 -1
//                             );
//
//                             const labelListList& targetFaces =
//                                 mpp.neighbPatch().AMI().tgtAddress();
//
//                             label subFaceI(0);
//                             forAll(pp.faceCells(), faceI)
//                             {
//                                 forAll(targetFaces[faceI], j)
//                                 {
//                                     faceBoundMap_[i][patchI][subFaceI] =
//                                         nbrNFaces;
//
//                                     subFaceI++;
//                                     nbrNFaces++;
//                                 }
//                             }
//                         }
//                     }
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
