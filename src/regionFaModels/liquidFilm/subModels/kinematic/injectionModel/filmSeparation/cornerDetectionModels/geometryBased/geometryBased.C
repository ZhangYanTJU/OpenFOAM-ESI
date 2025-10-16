/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "geometryBased.H"
#include "processorFaPatch.H"
#include "unitConversion.H"
#include "syncTools.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cornerDetectionModels
{
    defineTypeNameAndDebug(geometryBased, 0);
    addToRunTimeSelectionTable(cornerDetectionModel, geometryBased, dictionary);
}
}

const Foam::Enum
<
    Foam::cornerDetectionModels::geometryBased::cornerCurveType
>
Foam::cornerDetectionModels::geometryBased::cornerCurveTypeNames
({
    { cornerCurveType::ANY, "any" },
    { cornerCurveType::CONCAVE , "concave" },
    { cornerCurveType::CONVEX , "convex" }
});

const Foam::Enum
<
    Foam::cornerDetectionModels::geometryBased::cornerType
>
Foam::cornerDetectionModels::geometryBased::cornerTypeNames
({
    { cornerType::ALL, "sharpOrRound" },
    { cornerType::SHARP , "sharp" },
    { cornerType::ROUND , "round" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cornerDetectionModels::geometryBased::curvatureSign
(
    const vector& t,
    const vector& n0,
    const vector& n1
) const
{
    // t: unit edge tangent
    // n0: owner face unit normal
    // n1: neighbour face unit normal
    scalar curvature = (t & (n0 ^ n1));

    // Orientation: sign of triple product t . (n0 x n1)
    // Positive => one sense (together with outward normals, treat as "convex");
    // mapping to convex/concave is finally gated by 'cornerCurveType_'.
    return sign(curvature);
}


void Foam::cornerDetectionModels::geometryBased::classifyEdges
(
    bitSet& sharpEdges,
    bitSet& roundEdges
) const
{
    // Cache the operand references
    const areaVectorField& nf = mesh().faceAreaNormals();
    const edgeList& edges = mesh().edges();
    const labelUList& own = mesh().edgeOwner();      // own.sz = nEdges
    const labelUList& nei = mesh().edgeNeighbour();  // nei.sz = nInternalEdges
    const pointField& pts = mesh().points();


    // Convert input-angle parameters from degrees to radians
    const scalar angSharp = degToRad(angleSharpDeg_);
    const scalar angRoundMin = degToRad(angleRoundMinDeg_);
    const scalar angRoundMax = degToRad(angleRoundMaxDeg_);


    // Limit to subset of patches if requested
    bitSet allowedFaces(mesh().nFaces(), true);


    // Allocate the resource for the return objects
    sharpEdges.resize(mesh().nEdges());  // internal + boundary edges
    sharpEdges.reset();

    roundEdges.resize(mesh().nEdges());
    roundEdges.reset();


    // Internal edges
    const label nInternalEdges = mesh().nInternalEdges();
    for (label ei = 0; ei < nInternalEdges; ++ei)
    {
        // Do not allow processing of edges shorter than 'minEdgeLength'
        const edge& e = edges[ei];
        const scalar le = e.mag(pts);
        if (le <= max(minEdgeLength_, VSMALL)) continue;


        // Do not allow processing of excluded faces
        const label f0 = own[ei];
        const label f1 = nei[ei];
        if (!allowedFaces.test(f0) && !allowedFaces.test(f1)) continue;


        // Calculate the dihedral angle and curvature per edge
        const vector& n0 = nf[f0];
        const vector& n1 = nf[f1];

        const scalar phi = this->dihedralAngle(n0, n1); // [rad]
        const scalar kappa = 2.0*Foam::sin(0.5*phi)/max(le, VSMALL); // [1/m]
        const scalar R = (kappa > VSMALL ? scalar(1)/kappa : GREAT);

        const vector tangent(e.unitVec(pts));

        const scalar sgn = curvatureSign(tangent, n0, n1);

        const bool curvatureType =
            (cornerCurveType_ == cornerCurveType::ANY)
         || (cornerCurveType_ == cornerCurveType::CONVEX  && sgn > 0)
         || (cornerCurveType_ == cornerCurveType::CONCAVE && sgn < 0);

        // Do not allow processing of excluded curvature-type faces
        if (!curvatureType) continue;


        // Sharp: dihedral above threshold
        if (phi >= angSharp && cornerType_ != cornerType::ROUND)
        {
            sharpEdges.set(ei);
            continue; // do not double-classify as round
        }


        // Round: small-to-moderate angle but small radius (tight fillet)
        if
        (
            phi >= angRoundMin && phi <= angRoundMax
         && R <= maxRoundRadius_
         && cornerType_ != cornerType::SHARP
        )
        {
            roundEdges.set(ei);
        }
    }


    // Optional binary smoothing (edge-neighbour OR)


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return;

    // Boundary edges
    const faBoundaryMesh& patches = mesh().boundary();
    for (const faPatch& fap : patches)
    {
        const label patchi = fap.index();
        const label boundaryEdge0 = fap.start();

        const auto& nfp = nf.boundaryField()[patchi];

        if (isA<processorFaPatch>(fap))
        {
            forAll(nfp, bEdgei)
            {
                const label meshEdgei = boundaryEdge0 + bEdgei;

                // Do not allow processing of edges shorter than 'minEdgeLength'
                const edge& e = edges[meshEdgei];
                const scalar le = e.mag(pts);
                if (le <= max(minEdgeLength_, VSMALL)) continue;


                // Do not allow processing of excluded faces
                const label faceO = own[meshEdgei];
                if (!allowedFaces.test(faceO)) continue;


                // Fetch normal vector of owner and neigh faces
                const vector& n0 = nf[faceO];
                const vector& n1 = nfp[bEdgei];


                // Calculate the dihedral angle and curvature per edge
                const scalar phi = this->dihedralAngle(n0, n1); // [rad]
                const scalar kappa = 2.0*Foam::sin(0.5*phi)/max(le, VSMALL);
                const scalar R = (kappa > VSMALL ? scalar(1)/kappa : GREAT);

                const vector tangent(e.unitVec(pts));

                const scalar sgn = curvatureSign(tangent, n0, n1);

                const bool curvatureType =
                    (cornerCurveType_ == cornerCurveType::ANY)
                 || (cornerCurveType_ == cornerCurveType::CONVEX  && sgn > 0)
                 || (cornerCurveType_ == cornerCurveType::CONCAVE && sgn < 0);

                // Do not allow processing of excluded curvature-type faces
                if (!curvatureType) continue;


                // Sharp: dihedral above threshold
                if (phi >= angSharp && cornerType_ != cornerType::ROUND)
                {
                    sharpEdges.set(meshEdgei);
                    continue; // do not double-classify as round
                }


                // Round: small-to-moderate angle but small radius
                if
                (
                    phi >= angRoundMin && phi <= angRoundMax
                 && R <= maxRoundRadius_
                 && cornerType_ != cornerType::SHARP
                )
                {
                    roundEdges.set(meshEdgei);
                }
            }
        }
        else
        {
            forAll(nfp, bEdgei)
            {
                const label meshEdgei = boundaryEdge0 + bEdgei;
                const label faceO = own[meshEdgei];

                if (sharpBoundaryEdges_ && allowedFaces.test(faceO))
                {
                    sharpEdges.set(meshEdgei);
                }
                // Do not allow round edges on physical boundaries
            }
        }
    }
}


void Foam::cornerDetectionModels::geometryBased::edgesToFaces
(
    const bitSet& edgeMask,
    bitSet& faceMask
) const
{
    // Cache the operand references
    const labelUList& own = mesh().edgeOwner();
    const labelUList& nei = mesh().edgeNeighbour();


    // Allocate the resource for the return objects
    faceMask.resize(mesh().nFaces());
    faceMask.reset();


    // Internal edges
    const label nInternalEdges = mesh().nInternalEdges();
    for (label ei = 0; ei < nInternalEdges; ++ei)
    {
        if (edgeMask.test(ei))
        {
            // pick the intersecting owner and neighbour faces at the edge
            faceMask.set(nei[ei]);
        }
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return;


    // Boundary edges
    const faBoundaryMesh& patches = mesh().boundary();
    for (const faPatch& fap : patches)
    {
        const label bEdge0 = fap.start();
        const label nbEdges = fap.size();

        for (label bEdgei = 0; bEdgei < nbEdges; ++bEdgei)
        {
            const label meshEdgei = bEdge0 + bEdgei;

            if (edgeMask.test(meshEdgei))
            {
                faceMask.set(own[meshEdgei]);
            }
        }
    }
}


Foam::scalarList Foam::cornerDetectionModels::geometryBased::calcCornerAngles
(
    const bitSet& faces,
    const bitSet& edges
) const
{
    // Cache the operand references
    const areaVectorField& nf = mesh().faceAreaNormals();
    const labelUList& own = mesh().edgeOwner();
    const labelUList& nei = mesh().edgeNeighbour();


    // Allocate the resource for the return object
    scalarList cornerFaceAngles(mesh().faces().size(), Zero);


    // Internal edges
    const label nInternalEdges = mesh().nInternalEdges();
    for (label ei = 0; ei < nInternalEdges; ++ei)
    {
        if (!edges[ei]) continue;

        const label faceO = own[ei];
        const label faceN = nei[ei];

        // If neither adjacent face is flagged as a corner, skip the atan2 work
        if (!faces[faceO] && !faces[faceN]) continue;

        const scalar ang = this->dihedralAngle(nf[faceO], nf[faceN]);

        if (faces[faceO]) cornerFaceAngles[faceO] = ang;
        if (faces[faceN]) cornerFaceAngles[faceN] = ang;
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerFaceAngles;


    // Boundary edges
    const faBoundaryMesh& patches = mesh().boundary();
    for (const faPatch& fap : patches)
    {
        if (!isA<processorFaPatch>(fap)) continue;

        const label patchi = fap.index();
        const label bEdge0 = fap.start();

        const auto& nfp = nf.boundaryField()[patchi];

        forAll(nfp, bEdgei)
        {
            const label meshEdgei = bEdge0 + bEdgei;
            const label faceO = own[meshEdgei];

            // Only if the mesh edge is a corner and the owner face is a corner
            if (!edges[meshEdgei] || !faces[faceO]) continue;

            cornerFaceAngles[faceO] =
                this->dihedralAngle(nf[faceO], nfp[bEdgei]);
        }
    }

    return cornerFaceAngles;
}


void Foam::cornerDetectionModels::geometryBased::writeEdgesAndFaces() const
{
    // Cache the operand references
    const auto& edges = mesh().edges();
    const auto& faces = mesh().faces();
    const pointField& pts = mesh().points();

    // Generic writer for masked primitives (edge/face)
    auto writeMasked =
    [&](
        const auto& geom,
        const auto& mask,
        const char* file
    )
    {
        OBJstream os(mesh().time().path()/file);
        forAll(mask, i) if (mask[i]) os.write(geom[i], pts);
    };

    const bool writeSharp =
        (cornerType_ == cornerType::ALL || cornerType_ == cornerType::SHARP);
    const bool writeRound =
        (cornerType_ == cornerType::ALL || cornerType_ == cornerType::ROUND);

    if (writeSharp)
    {
        writeMasked(edges, sharpEdges_, "sharp-edges.obj");
        writeMasked(faces, sharpFaces_, "sharp-faces.obj");
    }

    if (writeRound)
    {
        writeMasked(edges, roundEdges_, "round-edges.obj");
        writeMasked(faces, roundFaces_, "round-faces.obj");
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cornerDetectionModels::geometryBased::geometryBased
(
    const faMesh& mesh,
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
    const dictionary& dict
)
:
    cornerDetectionModel(mesh, film, dict),
    init_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cornerDetectionModels::geometryBased::detectCorners()
{
    if (!init_ || mesh().moving())
    {
        // Identify and store sharp/round edges/faces
        classifyEdges(sharpEdges_, roundEdges_);
        edgesToFaces(sharpEdges_, sharpFaces_);
        edgesToFaces(roundEdges_, roundFaces_);


        // Collect all operand edges/faces
        cornerEdges_ = (sharpEdges_ | roundEdges_);
        cornerFaces_ = (sharpFaces_ | roundFaces_);


        // Pass the operand edges/faces to the film-separation model
        this->setCornerFaces(cornerFaces_);
        this->setCornerAngles
        (
            calcCornerAngles(cornerFaces_, cornerEdges_)
        );


        init_ = true;
    }


    // Write edges and faces as OBJ sets for debug purposes, if need be
    if (debug && mesh().time().writeTime())
    {
        writeEdgesAndFaces();
    }

    return true;
}


bool Foam::cornerDetectionModels::geometryBased::read(const dictionary& dict)
{
    if (!cornerDetectionModel::read(dict))
    {
        return false;
    }


    cornerCurveType_ = cornerCurveTypeNames.getOrDefault
    (
        "cornerCurveType",
        dict,
        cornerCurveType::ANY
    );

    cornerType_ = cornerTypeNames.getOrDefault
    (
        "cornerType",
        dict,
        cornerType::ALL
    );

    angleSharpDeg_ = dict.getOrDefault<scalar>("angleSharp", 45);
    angleRoundMinDeg_ = dict.getOrDefault<scalar>("angleRoundMin", 5);
    angleRoundMaxDeg_ = dict.getOrDefault<scalar>("angleRoundMax", 45);
    maxRoundRadius_ = dict.getOrDefault<scalar>("maxRoundRadius", 2e-3);

    minEdgeLength_ = dict.getOrDefault<scalar>("minEdgeLength", 0);
    nSmooth_ = dict.getOrDefault<label>("nSmooth", 0);

    sharpBoundaryEdges_ = dict.getOrDefault<bool>("sharpBoundaryEdges", false);


    // Validate the input parameters
    if (angleSharpDeg_ <= 0 || angleSharpDeg_ >= 180)
    {
        FatalIOErrorInFunction(dict)
            << "angleSharp (" << angleSharpDeg_
            << " deg) must be in (0, 180)."
            << exit(FatalIOError);
    }

    if
    (
        angleRoundMinDeg_ < 0 || angleRoundMaxDeg_ > 180
     || angleRoundMinDeg_ > angleRoundMaxDeg_
    )
    {
        FatalIOErrorInFunction(dict)
            << "Inconsistent round-angle range: angleRoundMin="
            << angleRoundMinDeg_ << " deg, angleRoundMax=" << angleRoundMaxDeg_
            << " deg. Require 0 <= min <= max <= 180."
            << exit(FatalIOError);
    }

    if (angleSharpDeg_ <= angleRoundMaxDeg_)
    {
        WarningInFunction
            << "angleSharp (" << angleSharpDeg_
            << " deg) <= angleRoundMax (" << angleRoundMaxDeg_
            << " deg): sharp vs round thresholds overlap; "
            << "classification may be ambiguous."
            << nl;
    }

    if (maxRoundRadius_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "maxRoundRadius must be non-negative."
            << exit(FatalIOError);
    }

    if (minEdgeLength_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "minEdgeLength must be non-negative."
            << exit(FatalIOError);
    }

    if (nSmooth_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "nSmooth must be non-negative."
            << exit(FatalIOError);
    }


    sharpEdges_.clear();
    roundEdges_.clear();
    cornerEdges_.clear();
    sharpFaces_.clear();
    roundFaces_.clear();
    cornerFaces_.clear();


    // Force the re-identification of corner edges/faces
    init_ = false;

    return true;
}

// ************************************************************************* //
