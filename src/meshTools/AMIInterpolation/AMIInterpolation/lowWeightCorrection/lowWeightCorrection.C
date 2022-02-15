/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "lowWeightCorrection.H"
#include "addToRunTimeSelectionTable.H"
#include "edgeTopoDistancesData.H"
#include "PatchEdgeFaceWave.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lowWeightCorrection, 0);
    addToRunTimeSelectionTable(AMIInterpolation, lowWeightCorrection, dict);
    addToRunTimeSelectionTable
    (
        AMIInterpolation,
        lowWeightCorrection,
        component
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lowWeightCorrection::findNearest
(
    const polyMesh& mesh,
    const primitivePatch& pp,
    const globalIndex& globalFaces,
    const labelList& uncoveredFaces,
    labelListList& nearestCoveredFaces
) const
{
    // Data on all edges and faces
    typedef edgeTopoDistancesData<label, primitivePatch> Type;

    List<Type> allEdgeInfo(pp.nEdges());
    List<Type> allFaceInfo(pp.size());


    // Seed all covered faces
    bitSet isUncovered(pp.size(), uncoveredFaces);

    DynamicList<label> changedEdges(pp.nEdges());
    DynamicList<Type> changedInfo(pp.nEdges());
    forAll(pp, facei)
    {
        if (!isUncovered(facei))
        {
            const labelList& fEdges = pp.faceEdges()[facei];

            const Type seed
            (
                labelList(1, 0),
                labelList(1, globalFaces.toGlobal(facei))
            );

            allFaceInfo[facei] = seed;
            for (const label edgei : fEdges)
            {
                changedEdges.append(edgei);
                changedInfo.append(seed);
            }
        }
    }


    // Walk
    Type::trackData td;
    td.n_ = nDonors_;

    PatchEdgeFaceWave
    <
        primitivePatch,
        Type,
        Type::trackData
    > calc
    (
        mesh,
        pp,
        changedEdges,
        changedInfo,
        allEdgeInfo,
        allFaceInfo,
        0,  //returnReduce(pp.nEdges(), sumOp<label>())
        td
    );
    calc.iterate(nIters_);  // should be enough iterations?


    nearestCoveredFaces.setSize(uncoveredFaces.size());
    forAll(uncoveredFaces, i)
    {
        const label facei = uncoveredFaces[i];
        if (allFaceInfo[facei].valid(calc.data()))
        {
            // Collect donor faces
            nearestCoveredFaces[i] = allFaceInfo[facei].data();
        }
    }
}


void Foam::lowWeightCorrection::calculateStencil
(
    const label facei,
    const point& sample,
    const labelList& covered,
    const labelListList& stencilAddressing,

    labelListList& addressing,
    scalarListList& weights,
    scalarField& sumWeights,
    const pointField& faceCentres
) const
{
    addressing[facei].clear();
    weights[facei].clear();
    sumWeights[facei] = Zero;

    // The (uncovered) facei has local donors

    // Add the donors
    for (const label coveredSlot : covered)
    {
        // Get remote slots for the local covered face
        const labelList& slots = stencilAddressing[coveredSlot];

        for (const label sloti : slots)
        {
            const point& donorPt = faceCentres[sloti];
            const label index = addressing[facei].find(sloti);
            if (index == -1)
            {
                addressing[facei].append(sloti);
                weights[facei].append(1.0/mag(sample-donorPt));
            }
            else
            {
                weights[facei][index] += 1.0/mag(sample-donorPt);
            }
        }
    }

    scalarList& wghts = weights[facei];
    const scalar w = sum(wghts);
    forAll(wghts, i)
    {
        wghts[i] /= w;
    }
    sumWeights[facei] = 1.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lowWeightCorrection::lowWeightCorrection
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    AMIInterpolation(dict, reverseTarget),
    lowWeightCorrection_(dict.get<scalar>("lowWeightCorrection")),
    nDonors_(dict.get<label>("nDonors")),
    nIters_(dict.get<label>("nIters")),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.subDict(type() + "Coeffs").get<word>("AMIMethod"),
            dict.subDict(type() + "Coeffs"),
            reverseTarget
        )
    )
{
    if (nDonors_ <= 0 || nIters_ <= 0)
    {
        WarningInFunction << "Disabled low-weight correction. Using "
            << AMIPtr_->type() << " AMI interpolation" << endl;
    }
}


Foam::lowWeightCorrection::lowWeightCorrection
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
:
    AMIInterpolation(requireMatch, reverseTarget, lowWeightCorrection),
    lowWeightCorrection_(-1),
    nDonors_(0),
    nIters_(0),
    AMIPtr_()
{}


Foam::lowWeightCorrection::lowWeightCorrection(const lowWeightCorrection& ami)
:
    AMIInterpolation(ami),
    lowWeightCorrection_(ami.lowWeightCorrection_),
    nDonors_(ami.nDonors_),
    nIters_(ami.nIters_),
    AMIPtr_(ami.AMIPtr_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lowWeightCorrection::calculate
(
    const polyMesh& mesh,
    const label srcPatchi,
    const primitivePatch& srcPatch,
    const label tgtPatchi,
    const primitivePatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (upToDate_)
    {
        return false;
    }


    AMIInterpolation::calculate(srcPatch, tgtPatch, surfPtr);
    AMIPtr_->calculate
    (
        mesh,
        srcPatchi,
        srcPatch,
        tgtPatchi,
        tgtPatch,
        surfPtr
    );

    // Take over AMI data
    // ~~~~~~~~~~~~~~~~~~

    upToDate_ = AMIPtr_->upToDate();
    // distributed / singlePatchProc
    singlePatchProc_ = AMIPtr_->singlePatchProc();;

    // Source patch
    srcMagSf_ = AMIPtr_->srcMagSf();
    srcAddress_ = AMIPtr_->srcAddress();
    srcWeights_ = AMIPtr_->srcWeights();
    srcWeightsSum_ = AMIPtr_->srcWeightsSum();
    srcCentroids_ = AMIPtr_->srcCentroids();
    //TBD: srcPatchPts_ = AMIPtr_->srcPatchPts();
    tsrcPatch0_ = AMIPtr_->srcPatch0();
    if (AMIPtr_->distributed())
    {
        srcMapPtr_.reset(new mapDistribute(AMIPtr_->srcMap()));
    }

    // Target
    tgtMagSf_ = AMIPtr_->tgtMagSf();
    tgtAddress_ = AMIPtr_->tgtAddress();
    tgtWeights_ = AMIPtr_->tgtWeights();
    tgtWeightsSum_ = AMIPtr_->tgtWeightsSum();
    tgtCentroids_ = AMIPtr_->tgtCentroids();
    //TBD: tgtPatchPts_ = AMIPtr_->tgtPatchPts();
    ttgtPatch0_ = AMIPtr_->tgtPatch0();
    if (AMIPtr_->distributed())
    {
        tgtMapPtr_.reset(new mapDistribute(AMIPtr_->tgtMap()));
    }


    // Walk out valid donors
    // ~~~~~~~~~~~~~~~~~~~~~

    if (nDonors_ > 0 && nIters_ > 0)
    {
        // Extend source side
        {
            // Low weight faces
            DynamicList<label> uncoveredFaces;
            forAll(srcWeightsSum(), facei)
            {
                if (srcWeightsSum()[facei] < lowWeightCorrection_)
                {
                    uncoveredFaces.append(facei);
                }
            }
            uncoveredFaces.shrink();


            // Global indexing for src faces
            const globalIndex globalFaces(srcPatch.size());


            // Find nearest face with high weight
            labelListList nearestCoveredFaces;
            findNearest
            (
                mesh,
                mesh.boundaryMesh()[srcPatchi],
                globalFaces,
                uncoveredFaces,
                nearestCoveredFaces
            );


            // Create map to get remote data into nearestCoveredFaces
            // order
            List<Map<label>> compactMap;
            const mapDistribute uncoveredToPatchMap
            (
                globalFaces,
                nearestCoveredFaces,
                compactMap
            );
            // Get srcPatch faceCentres in stencil ordering
            pointField stencilFcs(srcPatch.faceCentres());
            uncoveredToPatchMap.distribute(stencilFcs);
            // Get addressing (to donors) over in stencil ordering
            labelListList stencilAddressing(srcAddress());
            uncoveredToPatchMap.distribute(stencilAddressing);
            //scalarListList stencilWeights(srcWeights());
            //uncoveredToPatchMap.distribute(stencilWeights);


            // Target side patch centres
            tmp<pointField> totherFcs;
            if (AMIPtr_->distributed())
            {
                totherFcs = new pointField(tgtPatch.faceCentres());
                AMIPtr_->tgtMap().distribute(totherFcs.ref());
            }
            else
            {
                totherFcs = tgtPatch.faceCentres();
            }
            const pointField& otherFcs = totherFcs();


            //forAll(uncoveredFaces, i)
            //{
            //    const label facei = uncoveredFaces[i];
            //    const labelList& covered = nearestCoveredFaces[i];
            //    Pout<< "Uncovered face:" << facei
            //        << " low weihgt:" << srcWeightsSum()[facei]
            //        << " at:" << srcPatch.faceCentres()[facei]
            //        << " has local donors:" << endl;
            //    for (const label coveredSlot : covered)
            //    {
            //        Pout<< "    fc:" << stencilFcs[coveredSlot] << nl
            //            << "    which has remotes:"
            //            << stencilAddressing[coveredSlot] << nl
            //            << "    at:"
            //            <<  UIndirectList<point>
            //                (
            //                    otherFcs,
            //                    stencilAddressing[coveredSlot]
            //                ) << nl
            //            //<< "    with weights:"
            //            //<< stencilWeights[coveredSlot] << nl
            //            << endl;
            //    }
            //}


            // Re-do interpolation on uncovered faces
            forAll(uncoveredFaces, i)
            {
                const label facei = uncoveredFaces[i];

                calculateStencil
                (
                    facei,
                    srcPatch.faceCentres()[facei],
                    nearestCoveredFaces[i],
                    stencilAddressing,

                    srcAddress_,
                    srcWeights_,
                    srcWeightsSum_,
                    otherFcs
                );
            }

            //OBJstream os
            //(
            //    mesh.time().path()
            //   /mesh.boundaryMesh()[srcPatchi].name()+"_src.obj"
            //);
            //forAll(uncoveredFaces, i)
            //{
            //    const label facei = uncoveredFaces[i];
            //    const point& fc = srcPatch.faceCentres()[facei];
            //    const labelList& slots = srcAddress_[facei];
            //    for (const label sloti : slots)
            //    {
            //        os.write(linePointRef(fc, otherFcs[sloti]));
            //    }
            //}
        }


        // Extend target side
        {
            // Low weight faces
            DynamicList<label> uncoveredFaces;
            forAll(tgtWeightsSum(), facei)
            {
                if (tgtWeightsSum()[facei] < lowWeightCorrection_)
                {
                    uncoveredFaces.append(facei);
                }
            }
            uncoveredFaces.shrink();


            // Global indexing for tgt faces
            const globalIndex globalFaces(tgtPatch.size());


            // Find nearest face with high weight
            labelListList nearestCoveredFaces;
            findNearest
            (
                mesh,
                mesh.boundaryMesh()[tgtPatchi],
                globalFaces,
                uncoveredFaces,
                nearestCoveredFaces
            );


            // Create map to get remote data into nearestCoveredFaces
            // order
            List<Map<label>> compactMap;
            const mapDistribute uncoveredToPatchMap
            (
                globalFaces,
                nearestCoveredFaces,
                compactMap
            );
            // Get tgtPatch faceCentres in stencil ordering
            pointField stencilFcs(tgtPatch.faceCentres());
            uncoveredToPatchMap.distribute(stencilFcs);
            // Get addressing (to donors) over in stencil ordering
            labelListList stencilAddressing(tgtAddress());
            uncoveredToPatchMap.distribute(stencilAddressing);
            scalarListList stencilWeights(tgtWeights());
            uncoveredToPatchMap.distribute(stencilWeights);


            // Target side patch centres
            tmp<pointField> totherFcs;
            if (AMIPtr_->distributed())
            {
                totherFcs = new pointField(srcPatch.faceCentres());
                AMIPtr_->srcMap().distribute(totherFcs.ref());
            }
            else
            {
                totherFcs = srcPatch.faceCentres();
            }
            const pointField& otherFcs = totherFcs();


            //forAll(uncoveredFaces, i)
            //{
            //    const label facei = uncoveredFaces[i];
            //    const labelList& covered = nearestCoveredFaces[i];
            //    Pout<< "Uncovered face:" << facei
            //        << " low weihgt:" << tgtWeights()[facei]
            //        << " at:" << tgtPatch.faceCentres()[facei]
            //        << " has local donors:" << endl;
            //    for (const label coveredSlot : covered)
            //    {
            //        Pout<< "    fc:" << stencilFcs[coveredSlot] << nl
            //            << "    which has remotes:"
            //            << stencilAddressing[coveredSlot] << nl
            //            << "    at:"
            //            <<  UIndirectList<point>
            //                (
            //                    otherFcs,
            //                    stencilAddressing[coveredSlot]
            //                ) << nl
            //            << "    with weights:"
            //            << stencilWeights[coveredSlot] << nl
            //            << endl;
            //    }
            //}

            // Re-do interpolation on uncovered faces
            forAll(uncoveredFaces, i)
            {
                const label facei = uncoveredFaces[i];

                calculateStencil
                (
                    facei,
                    tgtPatch.faceCentres()[facei],
                    nearestCoveredFaces[i],
                    stencilAddressing,

                    tgtAddress_,
                    tgtWeights_,
                    tgtWeightsSum_,
                    otherFcs
                );
            }

            //OBJstream os
            //(
            //    mesh.time().path()
            //   /mesh.boundaryMesh()[srcPatchi].name()+"_tgt.obj"
            //);
            //forAll(uncoveredFaces, i)
            //{
            //    const label facei = uncoveredFaces[i];
            //    const point& fc = tgtPatch.faceCentres()[facei];
            //    const labelList& slots = tgtAddress_[facei];
            //    for (const label sloti : slots)
            //    {
            //        os.write(linePointRef(fc, otherFcs[sloti]));
            //    }
            //}
        }


        ////- Write face connectivity as OBJ file
        //writeFaceConnectivity
        //(
        //    srcPatch,
        //    tgtPatch,
        //    srcAddress_
        //);
    }

    return upToDate_;
}


void Foam::lowWeightCorrection::write(Ostream& os) const
{
    AMIInterpolation::write(os);
    os.writeEntry("nDonors", nDonors_);
    os.writeEntry("nIters", nIters_);
    os.beginBlock(word(this->type() + "Coeffs"));
    AMIPtr_->write(os);
    os.endBlock();
}


// ************************************************************************* //
