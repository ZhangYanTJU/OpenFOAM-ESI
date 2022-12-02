/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "demandDrivenData.H"
#include "faPatchFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(edgeInterpolation, 0);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::edgeInterpolation::clearOut() const
{
    lPNptr_.reset(nullptr);
    weightsPtr_.reset(nullptr);
    deltaCoeffsPtr_.reset(nullptr);
    nonOrthCorrectionVectorsPtr_.reset(nullptr);
    skewCorrectionVectorsPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::edgeInterpolation::edgeInterpolation(const faMesh& fam)
:
    faMesh_(fam),
    lPNptr_(nullptr),
    weightsPtr_(nullptr),
    deltaCoeffsPtr_(nullptr),
    nonOrthCorrectionVectorsPtr_(nullptr),
    skewCorrectionVectorsPtr_(nullptr),
    orthogonal_(false),
    skew_(true)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::edgeInterpolation::~edgeInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::edgeScalarField& Foam::edgeInterpolation::lPN() const
{
    if (!lPNptr_)
    {
        makeLPN();
    }

    return (*lPNptr_);
}


const Foam::edgeScalarField& Foam::edgeInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        makeWeights();
    }

    return (*weightsPtr_);
}


const Foam::edgeScalarField& Foam::edgeInterpolation::deltaCoeffs() const
{
    if (!deltaCoeffsPtr_)
    {
        makeDeltaCoeffs();
    }

    return (*deltaCoeffsPtr_);
}


bool Foam::edgeInterpolation::orthogonal() const
{
    if (orthogonal_ == false && !nonOrthCorrectionVectorsPtr_)
    {
        makeCorrectionVectors();
    }

    return orthogonal_;
}


const Foam::edgeVectorField&
Foam::edgeInterpolation::nonOrthCorrectionVector() const
{
    if (orthogonal())
    {
        FatalErrorInFunction
            << "cannot return correctionVectors; mesh is orthogonal"
            << abort(FatalError);
    }

    return (*nonOrthCorrectionVectorsPtr_);
}


bool Foam::edgeInterpolation::skew() const
{
    if (skew_ == true && !skewCorrectionVectorsPtr_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const Foam::edgeVectorField&
Foam::edgeInterpolation::skewCorrectionVectors() const
{
    if (!skew())
    {
        FatalErrorInFunction
            << "cannot return skewCorrectionVectors; mesh is now skewed"
            << abort(FatalError);
    }

    return (*skewCorrectionVectorsPtr_);
}


bool Foam::edgeInterpolation::movePoints() const
{
    orthogonal_ = false;
    skew_ = true;

    clearOut();

    return true;
}


const Foam::vector& Foam::edgeInterpolation::skewCorr(const label edgei) const
{
    #ifdef FA_SKEW_CORRECTION

    return
        (
            skewCorrectionVectorsPtr_
          ? (*skewCorrectionVectorsPtr_)[edgei]
          : pTraits<vector>::zero
        );

    #else

    return (*skewCorrectionVectorsPtr_)[edgei];

    #endif
}


void Foam::edgeInterpolation::makeLPN() const
{
    DebugInFunction
        << "Constructing geodesic distance between points P and N"
        << endl;

    const polyMesh& pMesh = mesh().mesh();

    lPNptr_.reset
    (
        new edgeScalarField
        (
            IOobject
            (
                "lPN",
                pMesh.pointsInstance(),
                pMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedScalar(dimLength, Zero)
        )
    );
    edgeScalarField& lPN = *lPNptr_;

    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    scalarField& lPNIn = lPN.primitiveFieldRef();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgei)
    {
        const vector& skewCorrEdge = skewCorr(edgei);

        scalar lPE =
            mag
            (
                edgeCentres[edgei]
              - skewCorrEdge
              - faceCentres[owner[edgei]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgei]]
              - edgeCentres[edgei]
              + skewCorrEdge
            );

        lPNIn[edgei] = (lPE + lEN);

        // Do not allow any mag(val) < SMALL
        if (mag(lPNIn[edgei]) < SMALL)
        {
            lPNIn[edgei] = SMALL;
        }
    }


    forAll(lPN.boundaryField(), patchi)
    {
        mesh().boundary()[patchi].makeDeltaCoeffs
        (
            lPN.boundaryFieldRef()[patchi]
        );

        lPN.boundaryFieldRef()[patchi] = 1.0/lPN.boundaryField()[patchi];
    }


    DebugInFunction
        << "Finished constructing geodesic distance PN"
        << endl;
}


void Foam::edgeInterpolation::makeWeights() const
{
    DebugInFunction
        << "Constructing weighting factors for edge interpolation"
        << endl;

    const polyMesh& pMesh = mesh().mesh();

    weightsPtr_.reset
    (
        new edgeScalarField
        (
            IOobject
            (
                "weights",
                pMesh.pointsInstance(),
                pMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedScalar(dimless, 1)
        )
    );
    edgeScalarField& weights = *weightsPtr_;


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    scalarField& weightsIn = weights.primitiveFieldRef();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgei)
    {
        const vector& skewCorrEdge = skewCorr(edgei);

        scalar lPE =
            mag
            (
                edgeCentres[edgei]
              - skewCorrEdge
              - faceCentres[owner[edgei]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgei]]
              - edgeCentres[edgei]
              + skewCorrEdge
            );

        // weight = (0,1]
        const scalar lPN = lPE + lEN;
        if (mag(lPN) > SMALL)
        {
            weightsIn[edgei] = lEN/lPN;
        }
    }

    forAll(mesh().boundary(), patchi)
    {
        mesh().boundary()[patchi].makeWeights
        (
            weights.boundaryFieldRef()[patchi]
        );
    }

    DebugInFunction
        << "Finished constructing weighting factors for face interpolation"
        << endl;
}


void Foam::edgeInterpolation::makeDeltaCoeffs() const
{
    DebugInFunction
        << "Constructing differencing factors array for edge gradient"
        << endl;

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    const polyMesh& pMesh = mesh().mesh();

    deltaCoeffsPtr_.reset
    (
        new edgeScalarField
        (
            IOobject
            (
                "deltaCoeffs",
                pMesh.pointsInstance(),
                pMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedScalar(dimless/dimLength, SMALL)
        )
    );
    edgeScalarField& deltaCoeffs = *deltaCoeffsPtr_;
    scalarField& dc = deltaCoeffs.primitiveFieldRef();

    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    const edgeVectorField& lengths = mesh().Le();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();

    // Calculate skewness correction vectors if necessary
    (void) skew();

    forAll(owner, edgei)
    {
        // Edge normal - area normal
        vector edgeNormal =
            normalised(lengths[edgei] ^ edges[edgei].vec(points));

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgei]]
          - faceCentres[owner[edgei]];

        unitDelta.removeCollinear(edgeNormal);
        unitDelta.normalise();


        const vector& skewCorrEdge = skewCorr(edgei);

        scalar lPE =
            mag
            (
                edgeCentres[edgei]
              - skewCorrEdge
              - faceCentres[owner[edgei]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgei]]
              - edgeCentres[edgei]
              + skewCorrEdge
            );

        scalar lPN = lPE + lEN;


        // Edge normal - area tangent
        edgeNormal = normalised(lengths[edgei]);

        // Do not allow any mag(val) < SMALL
        const scalar alpha = lPN*(unitDelta & edgeNormal);
        if (mag(alpha) > SMALL)
        {
            dc[edgei] = scalar(1)/max(alpha, 0.05*lPN);
        }
    }


    forAll(deltaCoeffs.boundaryField(), patchi)
    {
        mesh().boundary()[patchi].makeDeltaCoeffs
        (
            deltaCoeffs.boundaryFieldRef()[patchi]
        );
    }
}


void Foam::edgeInterpolation::makeCorrectionVectors() const
{
    DebugInFunction
        << "Constructing non-orthogonal correction vectors"
        << endl;

    const polyMesh& pMesh = mesh().mesh();

    nonOrthCorrectionVectorsPtr_.reset
    (
        new edgeVectorField
        (
            IOobject
            (
                "nonOrthCorrectionVectors",
                pMesh.pointsInstance(),
                pMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedVector(dimless, Zero)
        )
    );
    edgeVectorField& CorrVecs = *nonOrthCorrectionVectorsPtr_;

    // Set local references to mesh data
    const areaVectorField& faceCentres = mesh().areaCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const edgeVectorField& lengths = mesh().Le();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();

    scalarField deltaCoeffs(owner.size(), SMALL);

    vectorField& CorrVecsIn = CorrVecs.primitiveFieldRef();

    forAll(owner, edgei)
    {
        // Edge normal - area normal
        vector edgeNormal =
            normalised(lengths[edgei] ^ edges[edgei].vec(points));

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgei]]
          - faceCentres[owner[edgei]];

        unitDelta.removeCollinear(edgeNormal);
        unitDelta.normalise();

        // Edge normal - area tangent
        edgeNormal = normalised(lengths[edgei]);

        // Do not allow any mag(val) < SMALL
        const scalar alpha = unitDelta & edgeNormal;
        if (mag(alpha) > SMALL)
        {
            deltaCoeffs[edgei] = scalar(1)/alpha;
        }

        // Edge correction vector
        CorrVecsIn[edgei] =
            edgeNormal
          - deltaCoeffs[edgei]*unitDelta;
    }


    edgeVectorField::Boundary& CorrVecsbf = CorrVecs.boundaryFieldRef();

    forAll(CorrVecs.boundaryField(), patchi)
    {
        mesh().boundary()[patchi].makeCorrectionVectors(CorrVecsbf[patchi]);
    }


    constexpr scalar maxNonOrthoRatio = 0.1;
    scalar nonOrthoCoeff = 0;

    if (owner.size() > 0)
    {
        scalarField sinAlpha(deltaCoeffs*mag(CorrVecs.internalField()));

        forAll(sinAlpha, edgei)
        {
            sinAlpha[edgei] = max(-1, min(sinAlpha[edgei], 1));
        }

        nonOrthoCoeff = max(Foam::asin(sinAlpha)*radToDeg());
    }

    reduce(nonOrthoCoeff, maxOp<scalar>());

    DebugInFunction
        << "non-orthogonality coefficient = " << nonOrthoCoeff << " deg."
        << endl;

    if (nonOrthoCoeff < maxNonOrthoRatio)
    {
        nonOrthCorrectionVectorsPtr_.reset(nullptr);
    }

    orthogonal_ = !bool(nonOrthCorrectionVectorsPtr_);


    DebugInFunction
        << "Finished constructing non-orthogonal correction vectors"
        << endl;
}


void Foam::edgeInterpolation::makeSkewCorrectionVectors() const
{
    DebugInFunction
        << "Constructing skew correction vectors"
        << endl;

    const polyMesh& pMesh = mesh().mesh();

    skewCorrectionVectorsPtr_.reset
    (
        new edgeVectorField
        (
            IOobject
            (
                "skewCorrectionVectors",
                pMesh.pointsInstance(),
                pMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedVector(dimless, Zero)
        )
    );
    edgeVectorField& skewCorrVecs = *skewCorrectionVectorsPtr_;

    // Set local references to mesh data
    const areaVectorField& C = mesh().areaCentres();
    const edgeVectorField& Ce = mesh().edgeCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const pointField& points = mesh().points();
    const edgeList& edges = mesh().edges();


    forAll(neighbour, edgei)
    {
        const vector& P = C[owner[edgei]];
        const vector& N = C[neighbour[edgei]];
        const vector& S = points[edges[edgei].start()];

        // (T:Eq. 5.4)
        const vector d(N - P);
        const vector e(edges[edgei].vec(points));
        const vector de(d^e);
        const scalar alpha = magSqr(de);

        if (alpha < SMALL)
        {
            // Too small - skew correction remains zero
            continue;
        }
        const scalar beta = -((d^(S - P)) & de)/alpha;

        // (T:Eq. 5.3)
        const vector E(S + beta*e);

        skewCorrVecs[edgei] = Ce[edgei] - E;
    }


    edgeVectorField::Boundary& bSkewCorrVecs =
        skewCorrVecs.boundaryFieldRef();

    forAll(skewCorrVecs.boundaryField(), patchi)
    {
        faePatchVectorField& patchSkewCorrVecs = bSkewCorrVecs[patchi];

        if (patchSkewCorrVecs.coupled())
        {
            const labelUList& edgeFaces =
                mesh().boundary()[patchi].edgeFaces();

            const edgeList::subList patchEdges =
                mesh().boundary()[patchi].patchSlice(edges);

            vectorField ngbC(C.boundaryField()[patchi].patchNeighbourField());

            forAll(patchSkewCorrVecs, edgei)
            {
                const vector& P = C[edgeFaces[edgei]];
                const vector& N = ngbC[edgei];
                const vector& S = points[patchEdges[edgei].start()];

                // (T:Eq. 5.4)
                const vector d(N - P);
                const vector e(patchEdges[edgei].vec(points));
                const vector de(d^e);
                const scalar alpha = magSqr(de);

                if (alpha < SMALL)
                {
                    // Too small - skew correction remains zero
                    continue;
                }
                const scalar beta = -((d^(S - P)) & de)/alpha;

                const vector E(S + beta*e);

                patchSkewCorrVecs[edgei] =
                    Ce.boundaryField()[patchi][edgei] - E;
            }
        }
    }

    #ifdef FA_SKEW_CORRECTION

    constexpr scalar maxSkewRatio = 0.1;
    scalar skewCoeff = 0;

    forAll(own, edgei)
    {
        const scalar magSkew = mag(skewCorrVecs[edgei]);

        const scalar lPN =
            mag
            (
                Ce[edgei]
              - skewCorrVecs[edgei]
              - C[owner[edgei]]
            )
          + mag
            (
                C[neighbour[edgei]]
              - Ce[edgei]
              + skewCorrVecs[edgei]
            );

        const scalar ratio = magSkew/lPN;

        if (skewCoeff < ratio)
        {
            skewCoeff = ratio;

            if (skewCoeff > maxSkewRatio)
            {
                break;
            }
        }
    }

    DebugInFunction
        << "skew coefficient = " << skewCoeff << endl;

    if (skewCoeff < maxSkewRatio)
    {
        skewCorrectionVectorsPtr_.reset(nullptr);
    }

    #endif

    skew_ = bool(skewCorrectionVectorsPtr_);


    DebugInFunction
        << "Finished constructing skew correction vectors"
        << endl;
}


// ************************************************************************* //
