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

#include "fvMatrixAssembly.H"
#include "dictionary.H"
#include "GAMGSolver.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicFvPatchField.H"
#include "processorFvPatch.H"
#include "mappedWallPolyPatch.H"
#include "demandDrivenData.H"
#include "calculatedProcessorFvPatchField.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMatrixAssembly, 0);
}

//* * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrixAssembly::collectStencilData
(
    const tmpNrc<mapDistribute>& mapPtr,
    const labelListList& stencil,
    const Type& data,
    List<Type>& expandedData
)
{
    expandedData.setSize(stencil.size());

    if (mapPtr.valid())
    {
        Type work(data);
        mapPtr().distribute(work);

        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type::value_type>(work, slots)
            );
        }
    }
    else
    {
        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type::value_type>(data, slots)
            );
        }
    }
}


void Foam::fvMatrixAssembly::weights
(
    const List<pointField>& nbrCcs,
    const List<pointField>& nbrFaceCtr,
    const List<vectorField>& nbrNfs,
    const cyclicAMIFvPatch& pp,
    scalarField& w,
    vectorField& deltaCoeffs,
    const labelList& rmtFaces,
    const labelListList& sourceFaces
)
{
    scalarField nbrDeltas(w.size(), Zero);
    scalarField deltas(w.size(), Zero);

    vectorField remNbrPatchD(w.size(), Zero);
    vectorField remPatchD(w.size(), Zero);

    const pointField& faceCentres = pp.cyclicAMIPatch().faceCentres();
    const vectorField patchD(pp.coupledFvPatch::delta());
    const vectorField pnf(pp.nf());

    forAll (rmtFaces, subFaceI)
    {
        label faceI = rmtFaces[subFaceI];

        forAll(sourceFaces[faceI], j)
        {
            //nbr Remote
            nbrDeltas[subFaceI] =
                nbrNfs[faceI][j]
                & (nbrFaceCtr[faceI][j] - nbrCcs[faceI][j]);

            remNbrPatchD[subFaceI] =
                faceCentres[faceI] - nbrCcs[faceI][j];

            //local Remote
            deltas[subFaceI] = pnf[faceI] & patchD[faceI];

            remPatchD[subFaceI] = patchD[faceI];
        }
    }

    forAll(w, j)
    {
        scalar di = deltas[j];
        scalar dni = nbrDeltas[j];

        w[j] = dni/(di + dni + ROOTVSMALL);
        deltaCoeffs[j] = remPatchD[j] - remNbrPatchD[j];
    }

    Pout << "weights : " << w << endl;
    Pout << "deltaCoeffs : " << deltaCoeffs << endl;
}


void Foam::fvMatrixAssembly::addBoundaryDiag
(
    scalarField& diag
) const
{

    forAll(psis_, i)
    {
        //forAll(psi.boundaryField(), patchi)
        forAll(primitiveMesh_.patchMap()[i], patchi)
        {
            label globalPatchi = primitiveMesh_.patchMap()[i][patchi];

            if (globalPatchi != -1)
            {
                // Using any of the fvMatrix to addToInternal
                matrices_[0].addToInternalField
                (
                    primitiveMesh_.lduAddr().patchAddr(globalPatchi),
                    internalCoeffs_[globalPatchi],
                    diag
                );
            }
        }
    }
}


void Foam::fvMatrixAssembly::addBoundarySource
(
    Field<scalar>& source,
    const bool couples
) const
{
    forAll(psis_, i)
    {
        const volScalarField& psi = psis_[i];

        forAll(psi.boundaryField(), patchi)
        {
            const fvPatchField<scalar>& ptf = psi.boundaryField()[patchi];

            label globalPatchi = primitiveMesh_.patchMap()[i][patchi];

            // In global field mappedWall is not patch
            if (globalPatchi != -1)
            {
                const Field<scalar>& pbc = boundaryCoeffs_[globalPatchi];

                if (!ptf.coupled())
                {
                    matrices_[0].addToInternalField
                    (
                        primitiveMesh_.lduAddr().patchAddr(globalPatchi),
                        pbc,
                        source
                    );
                }
                else if (couples)
                {
                    const tmp<Field<scalar>> tpnf = ptf.patchNeighbourField();
                    const Field<scalar>& pnf = tpnf();

                    const labelUList& addr =
                        primitiveMesh_.lduAddr().patchAddr(globalPatchi);

                    forAll(addr, facei)
                    {
                        source[addr[facei]] +=
                            cmptMultiply(pbc[facei], pnf[facei]);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMatrixAssembly::fvMatrixAssembly
(
    const lduPrimitiveMeshAssembly& mesh,
    const dimensionSet& ds,
    const word& psiName
)
:
    lduMatrix(mesh),
    primitiveMesh_(mesh),
    matrices_(),
    psis_(),
    nMatrix_(0),
    psiName_(psiName),
    dimensions_(ds),
    source_(primitiveMesh_.lduAddr().size(), Zero),
    interfaces_(),
    faceAreasPtr_(nullptr),
    cellVolumesPtr_(nullptr)
{
    // Init upper and diagonal assuming symmetry matrix
    upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    diag().setSize(primitiveMesh_.lduAddr().size(), Zero);

    // We need the number of original patches + the new remoteInterfaces
    // To allocate internalCoeffs/boundaryCoeffs in order to map matrix.flux
    // from asembled to original
    label nPatches = primitiveMesh_.meshes()[0].boundary().size();
    for (label i=1; i < primitiveMesh_.meshes().size(); ++i)
    {
        nPatches += primitiveMesh_.meshes()[i].boundary().size();
    }
    nPatches += primitiveMesh_.remoteStencilInterfaces().size();

    // Size the internalCoeffs_/boundaryCoeffs_ with ALL the patches
    internalCoeffs_.setSize(nPatches);
    boundaryCoeffs_.setSize(nPatches);

    Info << "Total patches " << nPatches << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMatrixAssembly::calculateFlux()
{
    for (label i=0; i < nMatrix_; ++i)
    {
        // Use local flux() from the matrices and then modify the boundaries
        // flux
        tmp<surfaceScalarField> tflux = matrices_[i].flux();

        surfaceScalarField& flux = tflux.ref();

        // Fill internalCoeffs in flux
        forAll(matrices_[i].internalCoeffs(), patchi)
        {
            label globalPatchId = primitiveMesh_.patchMap()[i][patchi];

            if (globalPatchId == -1)
            {
                const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchi];

                if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);

                    if (mpp.owner())
                    {
                        const label nbrPatchId = mpp.neighbPatchID();

                        scalarField& internal =
                            matrices_[i].internalCoeffs()[patchi];

                        scalarField& nbrInternal =
                            matrices_[i].internalCoeffs()[nbrPatchId];

                        internal = Zero;
                        nbrInternal = Zero;

                        label virtualPatch =
                            primitiveMesh_.patchLocalToGlobalMap()[i][patchi];

                        const scalarField& internalAssembled =
                            internalCoeffs_[virtualPatch];

                        const labelList& cellIds =
                            primitiveMesh_.cellBoundMap()[i][nbrPatchId];

                        // nrb cells
                        const labelList& nbrCellIds =
                            primitiveMesh_.cellBoundMap()[i][patchi];

                        forAll(internalAssembled, subFaceI)
                        {
                            const label faceId =
                                primitiveMesh_.magSfFaceBoundMap()[i][patchi]
                                    [subFaceI];

                            const label nbrFaceId =
                                primitiveMesh_.magSfFaceBoundMap()[i][nbrPatchId]
                                [subFaceI];

                            if ((faceId != -1) && (nbrFaceId != -1))
                            {
                                const label cellId = cellIds[subFaceI];
                                const label nbrCellId = nbrCellIds[subFaceI];

                                internal[faceId] +=
                                    cmptMultiply
                                    (
                                        internalAssembled[subFaceI],
                                        psis_[i][cellId]
                                    );

                                nbrInternal[nbrFaceId] +=
                                    cmptMultiply
                                    (
                                        internalAssembled[subFaceI],
                                        psis_[i][nbrCellId]
                                    );
                            }
                        }
                    }
                }
            }
        }

        // Add contribution of new remote interfaces to original
        // cyclicAMI patches.

        const polyBoundaryMesh& patches = psis_[i].mesh().boundaryMesh();
        label firstRemote = patches.size();

        const PtrList<const lduPrimitiveProcessorInterface>& remoteInterfaces =
            primitiveMesh_.remoteStencilInterfaces();

        forAll(remoteInterfaces, remoteI)
        {
            label pint = primitiveMesh_.patchMap()[i][firstRemote++];

            const scalarField& assembledInternal = internalCoeffs_[pint];

            const labelUList& fc = remoteInterfaces[remoteI].faceCells();

            const scalarField internalField(psis_[i], fc);

            label supressed =
                primitiveMesh_.patchRemoteToLocal()[i][remoteI];

            // owner side
            const cyclicAMIPolyPatch& mpp =
                refCast<const cyclicAMIPolyPatch>(patches[supressed]);

            // check if remoteInterface & onwer in this proc has proc faces
            label neigProc = remoteInterfaces[remoteI].neighbProcNo();

            const Map<label>& map = mpp.AMI().tgtcMap()[neigProc];

            if (map.size() > 0)
            {
                // Do this side owner side on myProc that have proc faces on
                // the AMI patch
                scalarField& internal = matrices_[i].internalCoeffs()[supressed];

                forAll(assembledInternal, subFaceI)
                {
                    const label faceId =
                        primitiveMesh_.myRmtTgtFaces()
                        [i][mpp.index()][neigProc][subFaceI];

                    internal[faceId] +=
                        assembledInternal[subFaceI]*internalField[subFaceI];
                }
            }
            else
            {
                // Do the slave side on myProc
                scalarField& nbrInternal =
                    matrices_[i].internalCoeffs()[mpp.neighbPatchID()];

                forAll(assembledInternal, subFaceI)
                {
                    const label faceId =
                        primitiveMesh_.myRmtTgtFaces()[i][mpp.neighbPatchID()]
                        [neigProc][subFaceI];

                    nbrInternal[faceId] +=
                        assembledInternal[subFaceI]*internalField[subFaceI];
                }
            }
        }


        // Fill boundaryCoeffs
        forAll(matrices_[i].boundaryCoeffs(), patchi)
        {
            label globalPatchId = primitiveMesh_.patchMap()[i][patchi];

            if (globalPatchId == -1)
            {
                const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchi];

                if (isA<cyclicAMIPolyPatch>(pp))
                {
                     const cyclicAMIPolyPatch& mpp =
                        refCast<const cyclicAMIPolyPatch>(pp);

                    if (mpp.owner())
                    {
                        const label nbrPatchId = mpp.neighbPatchID();
                        scalarField& boundary = matrices_[i].boundaryCoeffs()
                            [patchi];
                        scalarField& nbrBoundary = matrices_[i].boundaryCoeffs()
                            [nbrPatchId];

                        boundary = Zero;
                        nbrBoundary = Zero;

                        label virtualPatch =
                            primitiveMesh_.patchLocalToGlobalMap()[i][patchi];

                        const scalarField& boundaryAssembled =
                            boundaryCoeffs_[virtualPatch];

                        // nrb cells
                        const labelList& nbrCellIds =
                            primitiveMesh_.cellBoundMap()[i][patchi];

                        const labelList& cellIds =
                            primitiveMesh_.cellBoundMap()[i][nbrPatchId];

                        forAll(boundaryAssembled, subFaceI)
                        {
                            const label faceId =
                                primitiveMesh_.magSfFaceBoundMap()[i][patchi]
                                    [subFaceI];

                            const label nbrFaceId =
                                primitiveMesh_.magSfFaceBoundMap()[i][nbrPatchId]
                                    [subFaceI];

                            if ((faceId != -1) && (nbrFaceId != -1))
                            {
                                const label nbrCellId = nbrCellIds[subFaceI];
                                const label cellId = cellIds[subFaceI];

                                boundary[faceId] +=
                                    cmptMultiply
                                    (
                                        boundaryAssembled[subFaceI],
                                        psis_[i][nbrCellId]
                                    );

                                nbrBoundary[nbrFaceId] +=
                                    cmptMultiply
                                    (
                                        boundaryAssembled[subFaceI],
                                        psis_[i][cellId]
                                    );
                            }
                        }
                    }
                }
            }
        }


        // Add rempte proc interfaces contribution to boundaryCoeffs
        firstRemote = patches.size();

        forAll(remoteInterfaces, remoteI)
        {
            const labelList& fc = remoteInterfaces[remoteI].faceCells();

            label pint = primitiveMesh_.patchMap()[i][firstRemote++];
            const scalarField& assembledBoundary = boundaryCoeffs_[pint];

            scalarField nbrInternalField(fc.size(), Zero);

            label supressed =
                primitiveMesh_.patchRemoteToLocal()[i][remoteI];

            const cyclicAMIPolyPatch& mpp =
                refCast<const cyclicAMIPolyPatch>(patches[supressed]);


            const labelListList& sourceFaces = mpp.AMI().srcAddress();
            const labelListList& tgtFaces = mpp.AMI().tgtAddress();
            const cyclicAMIPolyPatch& neigPatch = mpp.neighbPatch();

            const tmpNrc<mapDistribute> tgtMapPtr = mpp.AMI().tgtMap();

            // Get nbr internal values
            List<scalarField> srcFaceToTgtpsi;
            collectStencilData
            (
                tgtMapPtr,
                sourceFaces,
                scalarField(psis_[i], neigPatch.faceCells()),
                srcFaceToTgtpsi
            );
//DebugVar(srcFaceToTgtpsi)

            // Get my internal values
            List<scalarField> tgtFaceToSrcpsi;
            collectStencilData
            (
                mpp.AMI().srcMap(),
                tgtFaces,
                scalarField(psis_[i], mpp.faceCells()),
                tgtFaceToSrcpsi
            );


            // check if th9s remoteInterface & onwer in this proc has proc faces
            label neigProc = remoteInterfaces[remoteI].neighbProcNo();
            const labelList& srcMapSub = mpp.AMI().srcMap().subMap()[neigProc];

            label remoteFaceI = 0;
            forAll (srcMapSub, faceI)
            {
                const labelList& facesIds = sourceFaces[srcMapSub[faceI]];

                forAll(facesIds, j)
                {
                    label nbrFaceId = facesIds[j];
                    if (nbrFaceId >= neigPatch.size())
                    {
                        //Remote
                        nbrInternalField[remoteFaceI++] =
                            srcFaceToTgtpsi[srcMapSub[faceI]][j];
                    }
                }
            }
            if (srcMapSub.size() > 0)
            {
                scalarField& boundary = matrices_[i].boundaryCoeffs()[supressed];
                forAll(assembledBoundary, subFaceI)
                {
                    const label faceId =
                        primitiveMesh_.myRmtTgtFaces()
                            [i][mpp.index()][neigProc][subFaceI];

                    boundary[faceId] +=
                        assembledBoundary[subFaceI]*nbrInternalField[subFaceI];
                }
            }

            const labelList& tgtMapSub = mpp.AMI().tgtMap().subMap()[neigProc];
            remoteFaceI = 0;
            forAll (tgtMapSub, faceI)
            {
                const labelList& facesIds = tgtFaces[tgtMapSub[faceI]];

                forAll(facesIds, j)
                {
                    label faceId = facesIds[j];

                    if (faceId >= mpp.size())
                    {
                        nbrInternalField[remoteFaceI++] =
                                tgtFaceToSrcpsi[tgtMapSub[faceI]][j];
                    }
                }
            }

            if (tgtMapSub.size() > 0)
            {

                scalarField& nbrBoundary =
                    matrices_[i].boundaryCoeffs()[mpp.neighbPatchID()];

                forAll(assembledBoundary, subFaceI)
                {
                    const label faceId =
                        primitiveMesh_.myRmtTgtFaces()[i][mpp.neighbPatchID()]
                            [neigProc][subFaceI];

                    nbrBoundary[faceId] +=
                        assembledBoundary[subFaceI]*nbrInternalField[subFaceI];
                }
            }
        }

        surfaceScalarField::Boundary& ffbf = flux.boundaryFieldRef();

        forAll(ffbf, patchi)
        {
            label globalPatchId = primitiveMesh_.patchMap()[i][patchi];

            if (globalPatchId == -1)
            {
                ffbf[patchi] =
                    matrices_[i].internalCoeffs()[patchi]
                  - matrices_[i].boundaryCoeffs()[patchi];
            }

            //DebugVar(gSum(ffbf[patchi]))
        }

        matrices_[i].faceFluxPtr() = tflux.ptr();
    }

}


void Foam::fvMatrixAssembly::transferFieldsAndClean()
{
    // Reset local mesh
    //lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    //upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    //diag().setSize(primitiveMesh_.lduAddr().size(), Zero);
    //source_.setSize(primitiveMesh_.lduAddr().size(), Zero);

    const labelListList& procFaceMap = primitiveMesh_.faceMap();
    const labelList& cellMap = primitiveMesh_.cellOffsets();

    // Move append contents into intermediate list
    for (label i=0; i < nMatrix_; ++i)
    {
        if (asymmetric())
        {
            scalarField& lowerSub = matrices_[i].lower();
            forAll(lowerSub, facei)
            {
                lower()[procFaceMap[i][facei]] = lowerSub[facei];
            }
        }
        scalarField& upperSub = matrices_[i].upper();
        scalarField& diagSub = matrices_[i].diag();
        scalarField& sourceSub = matrices_[i].source();

        forAll(upperSub, facei)
        {
            upper()[procFaceMap[i][facei]] = upperSub[facei];
        }

        forAll(diagSub, celli)
        {
            const label globalCelli = cellMap[i] + celli;
            diag()[globalCelli] = diagSub[celli];
            source_[globalCelli] = sourceSub[celli];
        }

        if (!primitiveMesh_.fluxRequired(psiName_))
        {
            upperSub.clear();
            diagSub.clear();
            sourceSub.clear();
        }
    }
}


void Foam::fvMatrixAssembly::calcFaceAreasCellVolumes()
{
    const labelListList& faceMap = primitiveMesh_.faceMap();

    // Create the storage
    faceAreasPtr_ = new vectorField(lduAddr().upperAddr().size(), Zero);
    vectorField& faceAreas = *faceAreasPtr_;

    cellVolumesPtr_ = new scalarField(lduAddr().size(), Zero);
    scalarField& cellVolumes = *cellVolumesPtr_;

    for (label i=0; i < nMatrix_; ++i)
    {
        const scalarField& upperSub = matrices_[i].upper();
        const vectorField& areas = psis_[i].mesh().Sf();

        forAll(upperSub, facei)
        {
            faceAreas[faceMap[i][facei]] = areas[facei];
        }

        const polyBoundaryMesh& patches = psis_[i].mesh().boundaryMesh();

        // Fill faceAreas for new faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            label globalPatchID = primitiveMesh_.patchMap()[i][patchI];

            if (globalPatchID == -1)
            {
                if (isA<cyclicAMIPolyPatch>(pp))
                {
                    const cyclicAMIPolyPatch& mpp =
                            refCast<const cyclicAMIPolyPatch>(patches[patchI]);

                    if (mpp.owner())
                    {
                        const scalarListList& srcWeight =
                            mpp.AMI().srcWeights();

                        const vectorField& sf =
                            psis_[i].mesh().boundary()[patchI].Sf();


                        label subFaceI = 0;
                        forAll(pp.faceCells(), faceI)
                        {
                            const scalarList& w = srcWeight[faceI];

                            for(label j=0; j<w.size(); j++)
                            {
                                const label globalFaceI =
                                    primitiveMesh_.faceBoundMap()[i][patchI]
                                        [subFaceI];

                                if (globalFaceI != -1)
                                {
                                    faceAreas[globalFaceI] = w[j]*sf[faceI];
                                }
                                subFaceI++;
                            }

                        }
                    }
                }
            }
        }
        // Fill cellVolumes
        const scalarField& vol = psis_[i].mesh().V();
        const label cellOffset = primitiveMesh_.cellOffsets()[i];

        forAll(psis_[i], localCellI)
        {
            cellVolumes[cellOffset + localCellI] = vol[localCellI];
        }
    }
}


void Foam::fvMatrixAssembly::setInterfaces()
{
    if (interfaces_.empty())
    {
        interfaces_.setSize(internalCoeffs_.size());
        for (label i=0; i < nMatrix_; ++i)
        {
            const auto& bpsi = psis_[i].boundaryField();
            label interfaceID = 0;
            forAll(bpsi, patchI)
            {
                label globalPatchID =
                    primitiveMesh_.patchMap()[i][patchI];

                if (globalPatchID != -1)
                {
                    if (matrices_[i].mesh().interfaces().set(patchI))
                    {
                        if (isA<lduInterfaceField>(bpsi[patchI]))
                        {
                            if
                            (
                                isA<cyclicLduInterfaceField>(bpsi[patchI])
                             && primitiveMesh_.resetCyclics()
                            )
                            {
                                interfaces_.set
                                (
                                    interfaceID,
                                    new cyclicFvPatchField<scalar>
                                    (
                                        refCast<const fvPatch>
                                        (
                                            primitiveMesh_.interfaces()
                                            [
                                                globalPatchID
                                            ]
                                        ),
                                        bpsi[patchI].internalField()
                                    )
                                );
                            }
                            else
                            {
                                interfaces_.set
                                (
                                    interfaceID,
                                    &refCast<const lduInterfaceField>(bpsi[patchI])
                                );
                            }
                        }
                    }
                    interfaceID++;
                }
            }

            // Set remote interfaces
            forAll(primitiveMesh_.remoteStencilInterfaces(), remoteI)
            {

                label nOldInterfaces = matrices_[i].mesh().interfaces().size();
                label addPatchi = 0;
                for (label patchi = 0; patchi < nOldInterfaces; patchi++)
                {
                    if (isA<processorFvPatch>(bpsi[patchi].patch()))
                    {
                        addPatchi = patchi;
                        break;
                    }
                }

                interfaces_.set
                (
                    interfaceID,
                    new calculatedProcessorFvPatchField<scalar>
                    (
                        primitiveMesh_.interfaces()[interfaceID],
                        bpsi[addPatchi].patch(),    // dummy processorFvPatch
                        psis_[i]
                    )
                );
                interfaceID++;
            }
        }
    }
}

void Foam::fvMatrixAssembly::setBounAndInterCoeffs()
{
    for (label i=0; i < nMatrix_; ++i)
    {
        const polyBoundaryMesh& pbm = primitiveMesh_.meshes()[i].boundaryMesh();

        label interfaceID = 0;
        forAll(pbm, patchI)
        {
            label globalPatchId = primitiveMesh_.patchMap()[i][patchI];

            if (globalPatchId != -1)
            {
                if (matrices_[i].internalCoeffs().set(patchI))
                {
                    internalCoeffs_.set
                    (
                        globalPatchId,
                        matrices_[i].internalCoeffs()[patchI]
                    );
                }

                if (matrices_[i].boundaryCoeffs().set(patchI))
                {
                    boundaryCoeffs_.set
                    (
                        globalPatchId,
                        matrices_[i].boundaryCoeffs()[patchI]
                    );
                }
                interfaceID++;
            }
        }

        // Collect nbr data
        List<List<pointField>> nbrCcs(pbm.size());
        List<List<pointField>> nbrFaceCtrs(pbm.size());
        List<List<vectorField>> nbrNfs(pbm.size());

        if (primitiveMesh_.remoteStencilInterfaces().size() > 0)
        {
             // original interfaces
            const lduInterfacePtrsList meshInterfaces(psis_[i].mesh().interfaces());

            forAll(pbm, patchI)
            {
                label globalPatchId = primitiveMesh_.patchMap()[i][patchI];

                if (globalPatchId == -1)
                {
                    const polyPatch& pp = pbm[patchI];

                    if (isA<cyclicAMIPolyPatch>(pp))
                    {
                        const cyclicAMIPolyPatch& mpp =
                            refCast<const cyclicAMIPolyPatch>(pp);

                        if (mpp.owner())
                        {
                            const cyclicAMIFvPatch& pp =
                                refCast
                                <
                                    const cyclicAMIFvPatch
                                >(meshInterfaces[mpp.index()]);

                            const tmpNrc<mapDistribute> tgtMapPtr(mpp.AMI().tgtMap());

                            const cyclicAMIFvPatch& nbrPatch = pp.neighbPatch();

                            const labelListList& sourceFaces = mpp.AMI().srcAddress();

                            List<pointField> nbrCc;
                            collectStencilData
                            (
                                tgtMapPtr,
                                sourceFaces,
                                pointField(psis_[i].mesh().C(), nbrPatch.faceCells()),
                                nbrCc
                            );
                            nbrCcs[mpp.index()] = nbrCc;

                            List<pointField> nbrFaceCtr;
                            collectStencilData
                            (
                                tgtMapPtr,
                                sourceFaces,
                                pointField(nbrPatch.cyclicAMIPatch().faceCentres()),
                                nbrFaceCtr
                            );
                            nbrFaceCtrs[mpp.index()] = nbrFaceCtr;

                            List<vectorField> nbrNf;
                            const vectorField nf(nbrPatch.nf());
                            collectStencilData
                            (
                                tgtMapPtr,
                                sourceFaces,
                                nf,
                                nbrNf
                            );
                            nbrNfs[mpp.index()] = nbrNf;
                        }
                    }
                }
            }

        // Set new remote interfaces which follows

        // We need : gammaf*mag(Sf)*deltaCoeffs as boundaryCoeffs/internalCoeffs
        // for the new remote interfaces for the laplacian operator

        // weights: remote/local for gammaf interpolation
        // mag(Sf) on sub-face source local
        // deltaCoeffs : distance between cell centres

        // Normally all this is calculated by cyclicAMIFvPatch but here
        // we do it manually as fvPatch does not exist

        //if (primitiveMesh_.remoteStencilInterfaces().size() > 0)
        //{
            List<scalarField> internalCoeffs(Pstream::nProcs());
            List<scalarField> boundaryCoeffs(Pstream::nProcs());

            forAll(primitiveMesh_.remoteStencilInterfaces(), remoteI)
            {
                // original owner where remote originated
                label cyclicId =
                    primitiveMesh_.patchRemoteToLocal()[i][remoteI];

                const labelUList& fc =
                    primitiveMesh_.patchAddr()[interfaceID];

                const cyclicAMIFvPatch& pp =
                    refCast<const cyclicAMIFvPatch>(meshInterfaces[cyclicId]);

                label neigProc =
                    primitiveMesh_.remoteStencilInterfaces()
                    [
                        remoteI
                    ].neighbProcNo();

                const labelList& rmtFaces =
                    primitiveMesh_.myRmtTgtFaces()[0][pp.index()][neigProc];

                const labelListList& sourceFaces = pp.AMI().srcAddress();

                vectorField deltaCoeffs(fc.size(), Zero);
                scalarField w(fc.size(), Zero);

                // Calculate weights and deltaCoeffs
                weights
                (
                    nbrCcs[pp.index()],
                    nbrFaceCtrs[pp.index()],
                    nbrNfs[pp.index()],
                    pp,
                    w,
                    deltaCoeffs,
                    rmtFaces,
                    sourceFaces
                );

                internalCoeffs[Pstream::myProcNo()].setSize(fc.size(), Zero);
                boundaryCoeffs[Pstream::myProcNo()].setSize(fc.size(), Zero);

                internalCoeffs[neigProc].setSize(fc.size(), Zero);
                boundaryCoeffs[neigProc].setSize(fc.size(), Zero);

                // Calculate internalCoeffs/boundaryCoeffs on interfaces owner side
                // with connections to neigProc
                psis_[i].boundaryFieldRef()[cyclicId].manipulateInterBoundCoeffs
                (
                    *this,
                    w,
                    fc,
                    deltaCoeffs,
                    i,
                    neigProc,
                    internalCoeffs[Pstream::myProcNo()],
                    boundaryCoeffs[Pstream::myProcNo()]
                );

                Pstream::gatherList(internalCoeffs);
                Pstream::scatterList(internalCoeffs);

                Pstream::gatherList(boundaryCoeffs);
                Pstream::scatterList(boundaryCoeffs);

                const Map<label>& map = pp.AMI().tgtcMap()[neigProc];

                // From myProcNo to neigProc internalCoeffs > 0 for myProcNo
                // To make the internalCoeffs symmetrical copy to the other proc
                if (map.size() > 0)
                {
                    internalCoeffs_.set
                    (
                        interfaceID,
                        internalCoeffs[Pstream::myProcNo()]
                    );
                    boundaryCoeffs_.set
                    (
                        interfaceID,
                        boundaryCoeffs[Pstream::myProcNo()]
                    );
                }
                else
                {
                    internalCoeffs_.set
                    (
                        interfaceID,
                        internalCoeffs[neigProc]
                    );
                    boundaryCoeffs_.set
                    (
                        interfaceID,
                        boundaryCoeffs[neigProc]
                    );
                }

                Pout<< "internalCoeffs : interface : "
                    << interfaceID << " : " << internalCoeffs_[interfaceID]
                    << endl;

                interfaceID++;
            }
        }
    }
}

void Foam::fvMatrixAssembly::relaxFields()
{

    scalar relax(0.9);
    for (label i=0; i < nMatrix_; ++i)
    {
        forAll(psis_[i].mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchI];

            if (isA<cyclicAMIPolyPatch>(pp))
            {
                const labelList& cells = pp.faceCells();
                forAll(cells, facei)
                {
                    label celli = cells[facei];
                    psis_[i][celli] =
                        psis_[i][celli]*relax
                     +  psis_[i].oldTime()[celli]*(1 - relax);
                }
            }
        }
    }

}


void Foam::fvMatrixAssembly::manipulateMatrix()
{
    for (label i=0; i < nMatrix_; ++i)
    {
        forAll(psis_[i].mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchI];
            label globalPatchId = primitiveMesh_.patchMap()[i][patchI];

            if (globalPatchId == -1)
            {
                if (isA<mappedWallPolyPatch>(pp))
                {
                    volScalarField::Boundary& tbf =
                        matrices_[i].mesh().thisDb().lookupObjectRef
                        <
                            volScalarField
                        >("T").boundaryFieldRef();

                    tbf[patchI].manipulateMatrix
                    (
                        *this,
                        primitiveMesh_.faceBoundMap()[i][patchI],
                        primitiveMesh_.cellOffsets()[i],
                        i
                    );
                }
                else if (isA<cyclicAMIPolyPatch>(pp))
                {
                    psis_[i].boundaryFieldRef()[patchI].manipulateMatrix
                    (
                        *this,
                        primitiveMesh_.faceBoundMap()[i][patchI],
                        primitiveMesh_.cellOffsets()[i],
                        i
                    );
                }
                else
                {
                    FatalErrorInFunction
                        << "Local patch  : " << pp.index() << nl
                        << " it is not a mappedBased patch ID" << nl
                        << " in mesh : " << i
                        << exit(FatalError);
                }
            }
        }
    }
}


void Foam::fvMatrixAssembly::update()
{
    setInterfaces();
    // Transfer fields from matrices to local assembly
    // lower, upper, diag, source
    transferFieldsAndClean();
    setBounAndInterCoeffs();
    manipulateMatrix();
}


void Foam::fvMatrixAssembly::addFvMatrix(fvMatrix<scalar>& matrix)
{
    matrices_.append(&const_cast<fvMatrix<scalar>&>(matrix));
    psis_.append(&const_cast<volScalarField&>(matrix.psi()));

    // If one matrix is asymmetric the assembly matrix is asymmetric
    if (matrix.asymmetric())
    {
        lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    }
    ++nMatrix_;
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMatrixAssembly::solve
(
    const dictionary& solverControls
)
{
    const word solverType
    (
        solverControls.get<word>("solver")
    );

    if (solverType == GAMGSolver::typeName)
    {
        if (!faceAreasPtr_)
        {
            calcFaceAreasCellVolumes();
        }
    }

    // Update lower/upper/diag/source to this assembled matrix
    // and delete individual matrices
    update();

    scalarField saveDiag(diag());

    addBoundaryDiag(diag());

    addBoundarySource(source_, false);

    scalarField psi(lduAddr().size(), Zero);

    forAll(psis_, meshi)
    {
        label cellOffset = primitiveMesh_.cellOffsets()[meshi];

        forAll(psis_[meshi], localCellI)
        {
            psi[cellOffset + localCellI] = psis_[meshi][localCellI];
        }
    }

    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psiName_,
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces_,
        solverControls
    )->solve(psi, source_);


    forAll(psis_, i)
    {
        const label cellOffset = primitiveMesh_.cellOffsets()[i];

        forAll(psis_[i], localCellI)
        {
            psis_[i][localCellI] = psi[localCellI + cellOffset];
        }
    }


    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    // This will evaluate all BC's. The mapped BC should not alter the matrix.
    for (volScalarField& localPsi : psis_)
    {
        localPsi.correctBoundaryConditions();
        localPsi.mesh().setSolverPerformance(psiName_, solverPerf);
    }

    // Calculated locally the matrix flux
    if (primitiveMesh_.fluxRequired(psiName_))
    {
        calculateFlux();
    }
    //relaxFields();

    clear();

    return solverPerf;
}


void Foam::fvMatrixAssembly::clear()
{
    matrices_.clear();
    psis_.clear();
    if (asymmetric())
    {
        lower().clear();
    }
    upper().clear();
    diag().clear();
    nMatrix_ = 0;

    deleteDemandDrivenData(faceAreasPtr_);
    deleteDemandDrivenData(cellVolumesPtr_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMatrixAssembly::~fvMatrixAssembly()
{
    clear();
}


// ************************************************************************* //
