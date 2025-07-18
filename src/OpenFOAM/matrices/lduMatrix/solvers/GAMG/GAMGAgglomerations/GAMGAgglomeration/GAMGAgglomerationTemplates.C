/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023-2025 OpenCFD Ltd.
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

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const labelList& fineToCoarse
) const
{
    cf = Zero;

    forAll(ff, i)
    {
        cf[fineToCoarse[i]] += ff[i];
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    // By processor agglomeration the cf might be left zero size. Make sure
    // it is large enough.
    cf.resize_nocopy(nCells_[fineLevelIndex]);
    restrictField(cf, ff, fineToCoarse);

    const label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const auto& procIDs = agglomProcIDs(coarseLevelIndex);
        const auto& offsets = cellOffsets(coarseLevelIndex);

        globalIndex::gatherInplace
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = Zero;

    forAll(fineToCoarse, ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
            cf[cFace] += ff[ffacei];
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    const Field<Type>& cf,
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    const label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const auto& procIDs = agglomProcIDs(coarseLevelIndex);
        const auto& offsets = cellOffsets(coarseLevelIndex);

        const label localSize = nCells_[levelIndex];

        Field<Type> allCf(localSize);
        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            allCf,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );

        forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
}


template<class Type>
const Foam::Field<Type>& Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    Field<Type>& allCf,      // work storage
    const Field<Type>& cf,
    const label levelIndex
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    const label coarseLevelIndex = levelIndex+1;

    if (hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const auto& procIDs = agglomProcIDs(coarseLevelIndex);
        const auto& offsets = cellOffsets(coarseLevelIndex);

        const label localSize = nCells_[levelIndex];
        allCf.resize_nocopy(localSize);

        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            allCf,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );

        forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
        return allCf;
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
        return cf;
    }
}


// ************************************************************************* //
