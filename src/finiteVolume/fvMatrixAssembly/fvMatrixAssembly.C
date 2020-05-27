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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMatrixAssembly, 0);
}

//* * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::fvMatrixAssembly::weights
(
    const cyclicAMIFvPatch& pp,
    const labelUList& fc
)
{
    tmp<scalarField> tw(new scalarField(fc.size()));
    scalarField& w = tw.ref();

    const cyclicAMIPolyPatch& mpp = pp.cyclicAMIPatch();

    const cyclicAMIFvPatch& nbrPatch = pp.neighbFvPatch();

    const scalarField deltas(pp.nf() & pp.coupledFvPatch::delta());

    tmp<scalarField> tnbrDeltas;
        tnbrDeltas =
            mpp.interpolate
            (
                nbrPatch.nf()
              & nbrPatch.coupledFvPatch::delta()
            );

    const scalarField& nbrDeltas = tnbrDeltas();

    //primitiveMesh_.subFaceCompPatchMap()[i][mpp.index()]

    forAll(fc, j)
    {
        scalar di = deltas[fc[j]];
        scalar dni = nbrDeltas[fc[j]];

        w[j] = dni/(di + dni);
    }

    return tw;
}


Foam::tmp<Foam::vectorField> Foam::fvMatrixAssembly::deltaCoeffs
(
    const cyclicAMIFvPatch& pp,
    const labelUList& fc
)
{
    tmp<vectorField> tdeltaCoeffs(new vectorField(fc.size()));
    vectorField& deltaCoeffs = tdeltaCoeffs.ref();

    const cyclicAMIPolyPatch& mpp = pp.cyclicAMIPatch();

    const cyclicAMIFvPatch& nbrPatch = pp.neighbFvPatch();

    const vectorField patchD(pp.coupledFvPatch::delta());

    tmp<vectorField> tnbrPatchD;
    tnbrPatchD = mpp.interpolate(nbrPatch.coupledFvPatch::delta());
    const vectorField& nbrPatchD = tnbrPatchD();

    forAll(fc, j)
    {
        const vector& ddi = patchD[fc[j]];
        const vector& dni = nbrPatchD[fc[j]];

        deltaCoeffs[j] = ddi - dni;
    }

    return tdeltaCoeffs;
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

    label nPatches = primitiveMesh_.meshes()[0].boundary().size();
    for (label i=1; i < primitiveMesh_.meshes().size(); ++i)
    {
        nPatches += primitiveMesh_.meshes()[i].boundary().size();
    }
    forAll (primitiveMesh_.remoteStencilInterfaces(), remoteI)
    {
        if (primitiveMesh_.remoteStencilInterfaces().set(remoteI))
        {
            nPatches ++; //= primitiveMesh_.remoteStencilInterfaces().size();
        }
    }

    // Size the internalCoeffs_/boundaryCoeffs_ with ALL the patches
    internalCoeffs_.setSize(nPatches);
    boundaryCoeffs_.setSize(nPatches);

    DebugVar(nPatches);
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

                    label virtualPatch = -1;
                    if (mpp.owner())
                    {
                        const label nbrPatchId = mpp.neighbPatchID();
                        scalarField& internal = matrices_[i].internalCoeffs()
                            [patchi];
                        scalarField& nbrInternal = matrices_[i].internalCoeffs()
                            [nbrPatchId];

                        internal = Zero;
                        nbrInternal = Zero;

                        virtualPatch = primitiveMesh_.patchLocalToGlobalMap()[i]
                            [patchi];

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
//             else
//             {
//                 scalarField& internal = matrices_[i].internalCoeffs()[patchi];
//                 internal =
//                     cmptMultiply
//                     (
//                         internal,
//                         psis_[i].boundaryField()[patchi].patchInternalField()
//                     );
//             }
        }

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

                    label virtualPatch = -1;

                    if (mpp.owner())
                    {
                        const label nbrPatchId = mpp.neighbPatchID();
                        scalarField& boundary = matrices_[i].boundaryCoeffs()
                            [patchi];
                        scalarField& nbrBoundary = matrices_[i].boundaryCoeffs()
                            [nbrPatchId];

                        boundary = Zero;
                        nbrBoundary = Zero;

                        virtualPatch =
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

            DebugVar(gSum(ffbf[patchi]))
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
            //lowerSub.clear();
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

                        //const labelListList& sourceFaces =
                        //    mpp.AMI().srcAddress();

                        label subFaceI = 0;
                        forAll(pp.faceCells(), faceI)
                        {
                            const scalarList& w = srcWeight[faceI];
                            //const labelList& facesIds = sourceFaces[faceI];

                            for(label j=0; j<w.size(); j++)
                            {
                                //if (facesIds[j] < mpp.neighbPatch().size())
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
                else
                {

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

    //DebugVar(faceAreas)
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
            // Set interface if remote
DebugVar(interfaceID)
            forAll(primitiveMesh_.remoteStencilInterfaces(), remoteI)
            {
                if (primitiveMesh_.remoteStencilInterfaces().set(remoteI))
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

                    //label globalPatchID = primitiveMesh_.patchMap()[i][patchi];
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
DebugVar("Finish with setInterfaces")

}

void Foam::fvMatrixAssembly::setBounAndInterCoeffs()
{
    for (label i=0; i < nMatrix_; ++i)
    {
        label interfaceID = 0;
        forAll(psis_[i].mesh().boundaryMesh(), patchI)
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

        forAll(primitiveMesh_.remoteStencilInterfaces(), remoteI)
        {
            if (primitiveMesh_.remoteStencilInterfaces().set(remoteI))
            {
                label supressedPatch =
                    primitiveMesh_.patchRemoteToLocal()[i][remoteI];

                if (supressedPatch != -1)
                {
                    const labelUList& fc =
                        primitiveMesh_.patchAddr()[interfaceID];

// We need : gammaf*Sf*deltaCoeffs as boundaryCoeffs/internalCoeffs
// for the new processor Interface:
// weights: remote  for gamma interpolation on face
// Sf of source local
// deltaCoeffs : remote
// Normally all this is calculated by cyclicAMIFvPatch but here
// we do it manually as fvPatch does not exist
                    const cyclicAMIFvPatch& pp =
                        refCast<const cyclicAMIFvPatch>
                        (
                            psis_[i].boundaryField()[supressedPatch]
                        );

                    tmp<scalarField> w = weights(pp, fc);
                    //tmp<vectorField> localSf(pp.Sf(), fc);
                    tmp<vectorField> delta(deltaCoeffs(pp, fc));

                    scalarField internalCoeffs(fc.size(), Zero);
                    scalarField boundaryCoeffs(fc.size(), Zero);

                    psis_[i].boundaryFieldRef()
                    [
                        supressedPatch
                    ].manipulateInterBoundCoeffs
                    (
                        *this,
                        w,
                        fc,
                        delta,
                        i,
                        internalCoeffs,
                        boundaryCoeffs
                    );
                    internalCoeffs_.set
                    (
                        interfaceID,
                        internalCoeffs
                    );
                    boundaryCoeffs_.set
                    (
                        interfaceID,
                        boundaryCoeffs
                    );
                }
                interfaceID++;
            }
        }
    }
DebugVar("Finish with inter/bound coeffs")
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

    DebugVar("Finish with manipulate matrix")
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
DebugVar("update")
    scalarField saveDiag(diag());

    addBoundaryDiag(diag());
DebugVar("addBoundaryDiag")
    addBoundarySource(source_, false);
DebugVar("addBoundarySource")
    scalarField psi(lduAddr().size(), Zero);

    forAll(psis_, meshi)
    {
        label cellOffset = primitiveMesh_.cellOffsets()[meshi];

        forAll(psis_[meshi], localCellI)
        {
            psi[cellOffset + localCellI] = psis_[meshi][localCellI];
        }
    }
DebugVar("solver")
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
