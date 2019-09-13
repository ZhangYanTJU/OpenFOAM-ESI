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
#include "mappedWallPolyPatch.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMatrixAssembly, 0);
}

//* * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMatrixAssembly::addBoundaryDiag
(
    scalarField& diag
) const
{
    forAll(internalCoeffs_, patchi)
    {
        // Using any of the fvMatrix to addToInternal
        matrices_[0].addToInternalField
        (
            primitiveMesh_.lduAddr().patchAddr(patchi),
            internalCoeffs_[patchi],
            diag
        );
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
    internalCoeffs_(primitiveMesh_.patchAddr().size()),
    boundaryCoeffs_(internalCoeffs_.size()),
    interfaces_(),
    faceAreasPtr_(nullptr),
    cellVolumesPtr_(nullptr)
{
    lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    diag().setSize(primitiveMesh_.lduAddr().size(), Zero);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMatrixAssembly::transferFieldsAndClean()
{
    // Reset local mesh
    lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), Zero);
    diag().setSize(primitiveMesh_.lduAddr().size(), Zero);
    source_.setSize(primitiveMesh_.lduAddr().size(), Zero);

    const labelListList& procFaceMap = primitiveMesh_.faceMap();
    const labelList& cellMap = primitiveMesh_.cellOffsets();

    // Move append contents into intermediate list
    for (label i=0; i < nMatrix_; ++i)
    {
        scalarField& lowerSub = matrices_[i].lower();
        scalarField& upperSub = matrices_[i].upper();
        scalarField& diagSub = matrices_[i].diag();
        scalarField& sourceSub = matrices_[i].source();

        forAll(lowerSub, facei)
        {
            lower()[procFaceMap[i][facei]] = lowerSub[facei];
            upper()[procFaceMap[i][facei]] = upperSub[facei];
        }

        forAll(diagSub, celli)
        {
            const label globalCelli = cellMap[i] + celli;
            diag()[globalCelli] = diagSub[celli];
            source_[globalCelli] = sourceSub[celli];
        }

        lowerSub.clear();
        upperSub.clear();
        diagSub.clear();
        sourceSub.clear();
    }
}


void Foam::fvMatrixAssembly::calcFaceAreasCellVolumes()
{
    const labelListList& procFaceMap = primitiveMesh_.faceMap();

    // Create the storage
    faceAreasPtr_ = new vectorField(lduAddr().upperAddr().size());
    vectorField& faceAreas = *faceAreasPtr_;

    cellVolumesPtr_ = new scalarField(lduAddr().size());
    scalarField& cellVolumes = *cellVolumesPtr_;

    for (label i=0; i < nMatrix_; ++i)
    {
        scalarField& lowerSub = matrices_[i].lower();
        const vectorField& areas = psis_[i].mesh().Sf();

        forAll(lowerSub, facei)
        {
            faceAreas[procFaceMap[i][facei]] = areas[facei];
        }

        const polyBoundaryMesh& patches = psis_[i].mesh().boundaryMesh();

        // Fill faceAreas for new faces
        forAll(patches, patchI)
        {
            forAll(patches[patchI], faceI)
            {
                if (primitiveMesh_.faceBoundMap()[i][patchI].size() > 0)
                {
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>(patches[patchI]);

                    const label globalFaceI =
                        primitiveMesh_.faceBoundMap()[i][patchI][faceI];

                    if (globalFaceI != -1 && mpp.owner())
                    {
                        faceAreas[globalFaceI] =
                            psis_[i].mesh().boundary()[patchI].Sf()[faceI];
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


void Foam::fvMatrixAssembly::update()
{
    if (interfaces_.empty())
    {
        interfaces_.setSize(internalCoeffs_.size());
        for (label i=0; i < nMatrix_; ++i)
        {
            const auto& bpsi = psis_[i].boundaryField();

            forAll(bpsi, patchI)
            {
                if (matrices_[i].mesh().interfaces().set(patchI))
                {
                    if (isA<lduInterfaceField>(bpsi[patchI]))
                    {
                        label globalPatchID =
                            primitiveMesh_.patchMap()[i][patchI];
                        if (isA<cyclicLduInterfaceField>(bpsi[patchI]))
                        {
                            interfaces_.set
                            (
                                globalPatchID,
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
                                globalPatchID,
                                &refCast<const lduInterfaceField>(bpsi[patchI])
                            );
                        }
                    }
                }
            }
        }
    }

    // Transfer fields from matrices to local assembly
    // lower, upper, diag, source
    transferFieldsAndClean();

    for (label i=0; i < nMatrix_; ++i)
    {
        forAll(psis_[i].mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchI];
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

                    boundaryCoeffs_.set
                    (
                        globalPatchId,
                        matrices_[i].boundaryCoeffs()[patchI]
                    );

                }
            }
            else
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
                        primitiveMesh_.cellOffsets()[i]
                    );

                    //updateCoeffs(i, patchI);
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


void Foam::fvMatrixAssembly::addFvMatrix(const fvMatrix<scalar>& matrix)
{
    matrices_.append(matrix);
    psis_.append(&const_cast<volScalarField&>(matrix.psi()));
    ++nMatrix_;
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMatrixAssembly::solve
(
    const dictionary& solverControls
)
{
    const word solverType
    (
        solverControls.subDict(psiName_).get<word>("solver")
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
        solverControls.subDict(psiName_)
    )->solve(psi, source_);

    forAll(psis_, meshi)
    {
        const label cellOffset = primitiveMesh_.cellOffsets()[meshi];

        forAll(psis_[meshi], localCellI)
        {
            psis_[meshi][localCellI] = psi[localCellI + cellOffset];
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

    clear();

    return solverPerf;
}


void Foam::fvMatrixAssembly::clear()
{
    matrices_.clear();
    psis_.clear();
    lower().clear();
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
