/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "cyclicACMIGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        lduInterfaceField
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        cyclicACMIGAMGInterfaceField,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const cyclicAMILduInterfaceField& p =
        refCast<const cyclicAMILduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    Istream& is
)
:
    GAMGInterfaceField(GAMGCp, is),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(readBool(is)),
    rank_(readLabel(is))
{}


Foam::cyclicACMIGAMGInterfaceField::cyclicACMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& local,
    const UPtrList<lduInterfaceField>& other
)
:
    GAMGInterfaceField(GAMGCp, local),
    cyclicACMIInterface_(refCast<const cyclicACMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const auto& p = refCast<const cyclicACMILduInterfaceField>(local);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicACMIGAMGInterfaceField::ready() const
{
    if
    (
        UPstream::finishedRequests
        (
            recvRequests_.start(),
            recvRequests_.size()
        )
     && UPstream::finishedRequests
        (
            recvRequests1_.start(),
            recvRequests1_.size()
        )
    )
    {
        recvRequests_.clear();
        recvRequests1_.clear();

        if
        (
            UPstream::finishedRequests
            (
                sendRequests_.start(),
                sendRequests_.size()
            )
         && UPstream::finishedRequests
            (
                sendRequests1_.start(),
                sendRequests1_.size()
            )
        )
        {
            sendRequests_.clear();
            sendRequests1_.clear();
        }

        return true;
    }

    return false;
}


void Foam::cyclicACMIGAMGInterfaceField::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const auto& AMI =
    (
        cyclicACMIInterface_.owner()
      ? cyclicACMIInterface_.AMI()
      : cyclicACMIInterface_.neighbPatch().AMI()
    );

    if (AMI.distributed() && AMI.comm() != -1)
    {
        DebugPout<< "cyclicACMIFvPatchField::initInterfaceMatrixUpdate() :"
            << " interface:" << cyclicACMIInterface_.index()
            << " size:" << cyclicACMIInterface_.size()
            << " owner:" << cyclicACMIInterface_.owner()
            << " AMI distributed:" << AMI.distributed()
            << endl;

        // Start sending
        if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Can only evaluate distributed AMI with nonBlocking"
                << exit(FatalError);
        }

        // Get neighbouring field
        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIInterface_.neighbPatchID());

        solveScalarField pnf(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(pnf, cmpt);


        // Assert that all receives are known to have finished
        if (!recvRequests_.empty() || !recvRequests1_.empty())
        {
            FatalErrorInFunction
                << "Outstanding recv request(s) on patch "
                << cyclicACMIInterface_.index()
                << abort(FatalError);
        }

        const auto& cache = AMI.cache();

        if (cache.index0() == -1 && cache.index1() == -1)
        {
            const auto& map =
            (
                cyclicACMIInterface_.owner()
              ? AMI.tgtMap()
              : AMI.srcMap()
            );

            // Insert send/receive requests (non-blocking). See e.g.
            // cyclicAMIPolyPatchTemplates.C
            const label oldWarnComm = UPstream::commWarn(AMI.comm());
            map.send
            (
                pnf,
                sendRequests_,
                scalarSendBufs_,
                recvRequests_,
                scalarRecvBufs_,
                19462+cyclicACMIInterface_.index()  // unique offset + patch index
            );
            UPstream::commWarn(oldWarnComm);
        }
        else
        {
            cache.setDirection(cyclicACMIInterface_.owner());

            if (cache.index0() != -1)
            {
                const auto& map0 = cache.cTgtMapPtr0()();
                map0.send
                (
                    pnf,
                    sendRequests_,
                    scalarSendBufs_,
                    recvRequests_,
                    scalarRecvBufs_,
                    19462+cyclicACMIInterface_.index()   // unique offset + patch index
                );
            }

            if (cache.index1() != -1)
            {
                const auto& map1 = cache.cTgtMapPtr1()();
                map1.send
                (
                    pnf,
                    sendRequests1_,
                    scalarSendBufs1_,
                    recvRequests1_,
                    scalarRecvBufs1_,
                    19463+cyclicACMIInterface_.index()   // unique offset + patch index
                );
            }
        }
    }

    this->updatedMatrix(false);
}


void Foam::cyclicACMIGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    typedef multiplyWeightedOp<scalar, plusEqOp<scalar>> opType;

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    const auto& AMI =
    (
        cyclicACMIInterface_.owner()
      ? cyclicACMIInterface_.AMI()
      : cyclicACMIInterface_.neighbPatch().AMI()
    );

    const auto& cache = AMI.cache();

    if (AMI.distributed() && AMI.comm() != -1)
    {
        if (cache.index0() == -1 && cache.index1() == -1)
        {
            const auto& map =
            (
                cyclicACMIInterface_.owner()
              ? AMI.tgtMap()
              : AMI.srcMap()
            );

            // Receive (= copy) data from buffers into work. TBD: receive directly
            // into slices of work.
            solveScalarField work;
            map.receive
            (
                recvRequests_,
                scalarRecvBufs_,
                work,
                19462+cyclicACMIInterface_.index()  // unique offset + patch index
            );

            // Receive requests all handled by last function call
            recvRequests_.clear();

            solveScalarField pnf(faceCells.size(), Zero);
            AMI.weightedSum
            (
                cyclicACMIInterface_.owner(),
                work,
                pnf,               // result
                solveScalarField::null()
            );

            // Add result using coefficients
            this->addToInternalField(result, !add, faceCells, coeffs, pnf);
        }
        else
        {
            cache.setDirection(cyclicACMIInterface_.owner());

            solveScalarField pnf(faceCells.size());
            solveScalarField work(faceCells.size());

            if (cache.index0() != -1)
            {
                // Receive (= copy) data from buffers into work. TBD: receive directly
                // into slices of work.
                work = Zero;
                cache.cTgtMapPtr0()().receive
                (
                    recvRequests_,
                    scalarRecvBufs_,
                    work,
                    19462+cyclicACMIInterface_.index()   // unique offset + patch index
                );

                // Receive requests all handled by last function call
                recvRequests_.clear();

                pnf = Zero;
                AMIInterpolation::weightedSum
                (
                    AMI.lowWeightCorrection(),
                    cache.cSrcAddress0(),
                    cache.cSrcWeights0(),
                    cache.cSrcWeightsSum0(),
                    work,
                    opType(plusEqOp<scalar>()),
                    pnf,
                    solveScalarField::null()
                );

                pnf *= (1-cache.weight());

                // Add result using coefficients
                this->addToInternalField(result, !add, faceCells, coeffs, pnf);
            }

            if (cache.index1() != -1)
            {
                // Receive (= copy) data from buffers into work. TBD: receive directly
                // into slices of work.
                work = Zero;
                cache.cTgtMapPtr1()().receive
                (
                    recvRequests1_,
                    scalarRecvBufs1_,
                    work,
                    19463+cyclicACMIInterface_.index()   // unique offset + patch index
                );

                // Receive requests all handled by last function call
                recvRequests1_.clear();

                pnf = Zero;
                AMIInterpolation::weightedSum
                (
                    AMI.lowWeightCorrection(),
                    cache.cSrcAddress1(),
                    cache.cSrcWeights1(),
                    cache.cSrcWeightsSum1(),
                    work,
                    opType(plusEqOp<scalar>()),
                    pnf,
                    solveScalarField::null()
                );

                pnf *= cache.weight();

                // Add result using coefficients
                this->addToInternalField(result, !add, faceCells, coeffs, pnf);
            }
        }
    }
    else
    {
        // Get neighbouring field
        const labelUList& nbrFaceCells =
            lduAddr.patchAddr(cyclicACMIInterface_.neighbPatchID());

        solveScalarField work(psiInternal, nbrFaceCells);

        // Transform according to the transformation tensors
        transformCoupleField(work, cmpt);

        if (cache.index0() == -1 && cache.index1() == -1)
        {
            solveScalarField pnf(faceCells.size(), Zero);

            AMI.weightedSum
            (
                cyclicACMIInterface_.owner(),
                work,
                pnf,                // result
                solveScalarField::null()
            );

            const labelUList& faceCells = lduAddr.patchAddr(patchId);

            this->addToInternalField(result, !add, faceCells, coeffs, pnf);
        }
        else
        {
            cache.setDirection(cyclicACMIInterface_.owner());

            solveScalarField pnf(faceCells.size());

            if (cache.index0() != -1)
            {
                pnf = Zero;
                AMIInterpolation::weightedSum
                (
                    AMI.lowWeightCorrection(),
                    cache.cSrcAddress0(),
                    cache.cSrcWeights0(),
                    cache.cSrcWeightsSum0(),
                    work,
                    opType(plusEqOp<scalar>()),
                    pnf,
                    solveScalarField::null()
                );

                pnf *= (1 - cache.weight());

                this->addToInternalField(result, !add, faceCells, coeffs, pnf);
            }

            if (cache.index1() != -1)
            {
                pnf = Zero;
                AMIInterpolation::weightedSum
                (
                    AMI.lowWeightCorrection(),
                    cache.cSrcAddress1(),
                    cache.cSrcWeights1(),
                    cache.cSrcWeightsSum1(),
                    work,
                    opType(plusEqOp<scalar>()),
                    pnf,
                    solveScalarField::null()
                );

                pnf *= cache.weight();

                this->addToInternalField(result, !add, faceCells, coeffs, pnf);
            }
        }
    }

    this->updatedMatrix(true);
}


void Foam::cyclicACMIGAMGInterfaceField::write(Ostream& os) const
{
    //GAMGInterfaceField::write(os);
    os  << token::SPACE << doTransform()
        << token::SPACE << rank();
}


// ************************************************************************* //
