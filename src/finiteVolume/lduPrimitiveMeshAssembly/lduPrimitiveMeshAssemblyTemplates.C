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

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::lduPrimitiveMeshAssembly::assemble()
{
    // Get newFaces and newPatches (not mapPolyPatches)
    label oldFaces(0);
    label newFaces(0);
    label newPatches(0);
    label nMeshes = meshes_.size();

    for (label i=0; i < nMeshes; ++i)
    {
        patchMap_[i].setSize(meshes_[i].boundaryMesh().size(), -1);

        patchLocalToGlobalMap_[i].setSize(meshes_[i].boundaryMesh().size(), -1);

        forAll(meshes_[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes_[i].boundaryMesh()[patchI];

            if (!isA<Type>(pp))
            {
                patchMap_[i][patchI] = newPatches;
                patchLocalToGlobalMap_[i][patchI] = newPatches;
                newPatches++;
            }
            else
            {
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
                        const scalarListList& weightSourceFaces =
                            mpp.AMI().srcWeights();

                        //newFaces += pp.size();

                        forAll (weightSourceFaces, faceI)
                        {
                            newFaces += weightSourceFaces[faceI].size();
                        }
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "The patch" << pp.name() << " is not of type "
                        << Type::typeName << " But is not of type "
                        << mappedPatchBase::typeName << " or "
                        << cyclicAMIPolyPatch::typeName
                        << " which are the two types which can be interlalized "
                        << exit(FatalError);
                }
            }
        }
    }

    label virtualPatches = newPatches;
    DebugVar(patchLocalToGlobalMap_)
    DebugVar(patchMap_)
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

    Info<< "old nFaces : " << oldFaces
        << "new nFaces : " << newFaces << endl;

    // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    Info<< "new patches : " << newPatches << endl;

    //DebugVar(patchMap_)

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
    }

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

    // Add lower/upper adressing for new internal faces corresponding
    // to old mapPolyPatch's
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
                    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

                    if (mpp.owner())
                    {
                        cellBoundMap_[i][patchI].setSize(pp.size(), -1);

                        nbrFaceCells =
                        (
                            mpp.samplePolyPatch().faceCells()
                        );

                        const label nbrPatchId = mpp.samplePolyPatch().index();
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
                        const labelListList& sourceFaces = mpp.AMI().srcAddress();

                        nbrFaceCells =
                        (
                            mpp.neighbPatch().faceCells()
                        );

                        const label nbrPatchId = mpp.neighbPatchID();

                        label newFaces(0);
                        forAll (sourceFaces, faceI)
                        {
                            newFaces += sourceFaces[faceI].size();
                        }


                        cellBoundMap_[i][patchI].setSize(newFaces, -1);
                        cellBoundMap_[i][nbrPatchId].setSize(newFaces, -1);
                        magSfFaceBoundMap_[i][patchI].setSize(newFaces, -1);
                        magSfFaceBoundMap_[i][nbrPatchId].setSize(newFaces, -1);

                        label subFaceI(0);
                        forAll(pp.faceCells(), faceI)
                        {

                            const label cellI =
                                pp.faceCells()[faceI] + cellOffsets_[i];

                            //const scalarList& wSrcFaces =  weightSourceFaces[faceI];
                            //const scalarList& wTgtFaces = weightTargetFaces[faceI];

                            const labelList& facesIds = sourceFaces[faceI];

                            forAll(facesIds, j)
                            {
                                label nbrFaceId = facesIds[j];

                                const label nbrCellI =
                                    nbrFaceCells[nbrFaceId] + cellOffsets_[meshNrbId];

                                lowerAddr()[nFaces] = min(cellI, nbrCellI);
                                upperAddr()[nFaces] = max(cellI, nbrCellI);

                                cellBoundMap_[i][patchI][subFaceI] = nbrCellI;
                                cellBoundMap_[i][nbrPatchId][subFaceI] = cellI;

                                magSfFaceBoundMap_[i][patchI][subFaceI] = faceI;
                                //magSfFaceBoundMap_[i][nbrPatchId][subFaceI] = faceI;

                                ++subFaceI;
                                ++nFaces;
                            }
                        }
                     }
                     else
                     {
                        const labelListList& tgtFaces =
                            mpp.neighbPatch().AMI().tgtAddress();
                        label subFaceI(0);

                        forAll(pp.faceCells(), faceI)
                        {
                            const labelList& facesIds = tgtFaces[faceI];
                            forAll(facesIds, j)
                            {
                                magSfFaceBoundMap_[i][patchI][subFaceI] = faceI;
                                ++subFaceI;
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
            if (isA<Type>(pp))// && !pp.empty())
            {
                label meshNrbId = 0;
                label samplePatchi = -1;
                //label ppNbrSize = 0;

                if (isA<mappedPatchBase>(pp))
                {
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>(pp);
                    if (mpp.owner())
                    {
                        samplePatchi = mpp.samplePolyPatch().index();
                        meshNrbId = this->findNbrMeshId(mpp, meshes_);
                        //const polyPatch& ppNbr = mpp.samplePolyPatch();
                        //ppNbrSize = ppNbr.size();

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
                            faceBoundMap_[meshNrbId][samplePatchi][faceI] = nFaces;
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

                        if
                        (
                            faceBoundMap_[i][patchI].empty()
                         && faceBoundMap_[i][samplePatchi].empty()
                        )
                        {
                            faceBoundMap_[i][patchI].setSize
                            (
                                cellBoundMap_[i][patchI].size(),
                                -1
                            );
                            faceBoundMap_[i][samplePatchi].setSize
                            (
                                cellBoundMap_[i][patchI].size(),
                                -1
                            );

                            const labelListList& sourceFaces = mpp.AMI().srcAddress();

                            label subFaceI(0);
                            forAll(pp.faceCells(), faceI)
                            {
                                forAll(sourceFaces[faceI], j)
                                {
                                    faceBoundMap_[i][patchI][subFaceI] = nFaces;
                                    faceBoundMap_[i][samplePatchi][subFaceI] = nFaces;
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

    if (newFaces != nFaces)
    {
       FatalErrorInFunction
            << "Incorrrect total number of faces in the assembled lduMatrix: "
            << newFaces << " != " << nFaces << nl
            << exit(FatalError);
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

//         for (labelListList& bMap : magSfFaceBoundMap_)
//         {
//             for (labelList& faceMap : bMap)
//             {
//                 for (label& facei : faceMap)
//                 {
//                     facei = oldToNew[facei];
//                 }
//             }
//         }

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
