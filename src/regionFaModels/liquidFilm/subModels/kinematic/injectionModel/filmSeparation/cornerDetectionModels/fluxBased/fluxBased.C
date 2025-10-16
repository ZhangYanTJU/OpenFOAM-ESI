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

#include "fluxBased.H"
#include "OBJstream.H"
#include "processorFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cornerDetectionModels
{
    defineTypeNameAndDebug(fluxBased, 0);
    addToRunTimeSelectionTable(cornerDetectionModel, fluxBased, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::bitSet Foam::cornerDetectionModels::fluxBased::identifyCornerEdges() const
{
    // Return true if face normals converge, i.e. sharp edge
    // Face-normal vectors diverge: no separation, converge: separation (maybe)
    const auto isCornerEdgeSharp =
    [](
        const vector& fcO,  // face-centre owner
        const vector& fcN,  // face-centre neigh
        const vector& fnO,  // face-normal owner
        const vector& fnN   // face-normal neigh
    ) noexcept -> bool
    {
        // Threshold for sharpness detection
        constexpr scalar sharpEdgeThreshold = -1e-8;

        // Relative centre and normal of the two faces sharing the edge
        const vector relativePosition(fcN - fcO);
        const vector relativeNormal(fnN - fnO);

        // Sharp if normals converge along the centre-to-centre direction
        return ((relativeNormal & relativePosition) < sharpEdgeThreshold);
    };


    // Cache the operand references
    const areaVectorField& fc = mesh().areaCentres();
    const areaVectorField& fn = mesh().faceAreaNormals();
    const labelUList& own = mesh().edgeOwner();
    const labelUList& nei = mesh().edgeNeighbour();


    // Allocate the resource for the return object
    bitSet cornerEdges(mesh().nEdges(), false);

    // Internal edges (owner <-> neighbour)
    forAll(nei, edgei)
    {
        const label faceO = own[edgei];
        const label faceN = nei[edgei];

        cornerEdges[edgei] = isCornerEdgeSharp
        (
            fc[faceO],
            fc[faceN],
            fn[faceO],
            fn[faceN]
        );
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerEdges;

    // Check if processor face-normal vectors diverge (no separation)
    // or converge (separation may occur)
    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (!isA<processorFaPatch>(fap)) continue;

        const label patchi = fap.index();
        const auto& edgeFaces = fap.edgeFaces();
        const label internalEdge0 = fap.start();

        const auto& fcp = fc.boundaryField()[patchi];
        const auto& fnp = fn.boundaryField()[patchi];

        // Processor edges (owner <-| none)
        forAll(fnp, bndEdgei)
        {
            const label faceO = edgeFaces[bndEdgei];
            const label meshEdgei = internalEdge0 + bndEdgei;

            cornerEdges[meshEdgei] = isCornerEdgeSharp
            (
                fc[faceO],
                fcp[bndEdgei],
                fn[faceO],
                fnp[bndEdgei]
            );
        }
    }

    return cornerEdges;
}


Foam::bitSet Foam::cornerDetectionModels::fluxBased::identifyCornerFaces
(
    const bitSet& cornerEdges
) const
{
    // Marks the separating face based on edge flux sign
    const auto markSeparation =
    [](
        bitSet& cornerFaces,
        const scalar phiEdge,
        const label faceO,
        const label faceN = -1  /* = -1 for processor edges */
    ) noexcept -> void
    {
        constexpr scalar tol = 1e-8;

        // Assuming no sources/sinks at the edge
        if (phiEdge > tol)  // From owner to neighbour
        {
            cornerFaces[faceO] = true;
        }
        else if ((phiEdge < -tol) && (faceN != -1))  // From nei to own
        {
            cornerFaces[faceN] = true;
        }
    };


    // Cache the operand references
    const edgeScalarField& phis = film().phi2s();
    const labelUList& own = mesh().edgeOwner();
    const labelUList& nei = mesh().edgeNeighbour();

    // Allocate the resource for the return object
    bitSet cornerFaces(mesh().faces().size(), false);

    // Internal faces (owner <-> neighbour)
    forAll(nei, edgei)
    {
        if (!cornerEdges[edgei]) continue;

        markSeparation
        (
            cornerFaces,
            phis[edgei],
            own[edgei],  // faceO
            nei[edgei]   // faceN
        );
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerFaces;

    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (!isA<processorFaPatch>(fap)) continue;

        const label patchi = fap.index();
        const auto& edgeFaces = fap.edgeFaces();
        const label internalEdge0 = fap.start();

        const auto& phisp = phis.boundaryField()[patchi];

        // Processor faces (owner <-| none)
        forAll(phisp, bndEdgei)
        {
            const label faceO = edgeFaces[bndEdgei];
            const label meshEdgei = internalEdge0 + bndEdgei;

            if (!cornerEdges[meshEdgei]) continue;

            markSeparation
            (
                cornerFaces,
                phisp[bndEdgei],
                faceO
                /*faceN = -1*/
            );
        }
    }

    return cornerFaces;
}


Foam::scalarList Foam::cornerDetectionModels::fluxBased::calcCornerAngles
(
    const bitSet& faces,
    const bitSet& edges
) const
{
    // Cache the operand references
    const areaVectorField& fn = mesh().faceAreaNormals();
    const labelUList& own = mesh().edgeOwner();
    const labelUList& nei = mesh().edgeNeighbour();

    scalarList cornerFaceAngles(mesh().faces().size(), Zero);

    // Internal edges (owner <-> neighbour)
    forAll(nei, edgei)
    {
        if (!edges[edgei]) continue;

        const label faceO = own[edgei];
        const label faceN = nei[edgei];

        // If neither adjacent face is flagged as a corner, skip the atan2 work
        if (!faces[faceO] && !faces[faceN]) continue;

        const scalar ang = this->dihedralAngle(fn[faceO], fn[faceN]);

        if (faces[faceO]) cornerFaceAngles[faceO] = ang;
        if (faces[faceN]) cornerFaceAngles[faceN] = ang;
    }


    // Skip the rest of the routine if the simulation is a serial run
    if (!Pstream::parRun()) return cornerFaceAngles;

    const faBoundaryMesh& patches = mesh().boundary();

    for (const faPatch& fap : patches)
    {
        if (!isA<processorFaPatch>(fap)) continue;
        const label patchi = fap.index();
        const auto& edgeFaces = fap.edgeFaces();
        const label internalEdge0 = fap.start();

        const auto& fnp = fn.boundaryField()[patchi];

        // Processor edges (owner <-| none)
        forAll(fnp, bndEdgei)
        {
            const label faceO = edgeFaces[bndEdgei];
            const label meshEdgei = internalEdge0 + bndEdgei;

                // Only if the mesh edge and owner face are both corners
                if (!edges[meshEdgei] || !faces[faceO]) continue;

                cornerFaceAngles[faceO] =
                    this->dihedralAngle(fn[faceO], fnp[bndEdgei]);
        }
    }

    return cornerFaceAngles;
}


void Foam::cornerDetectionModels::fluxBased::writeEdgesAndFaces
(
    const word& prefix
) const
{
    const pointField& pts = mesh().points();

    const word timeName(Foam::name(mesh().time().value()));
    const word nameEdges("fluxBased-edges-" + timeName + ".obj");
    const word nameFaces("fluxBased-faces-" + timeName + ".obj");


    // Write OBJ of edge faces to file
    OBJstream osEdges(mesh().time().path()/nameEdges);

    const auto& edges = mesh().edges();
    forAll(cornerEdges_, ei)
    {
        if (cornerEdges_[ei])
        {
            const edge& e = edges[ei];
            osEdges.write(e, pts);
        }
    }

    // Write OBJ of corner faces to file
    OBJstream osFaces(mesh().time().path()/nameFaces);

    const bitSet& cornerFaces = this->getCornerFaces();
    const auto& faces = mesh().faces();
    forAll(cornerFaces, fi)
    {
        if (cornerFaces[fi])
        {
            const face& f = faces[fi];
            osFaces.write(f, pts);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cornerDetectionModels::fluxBased::fluxBased
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

bool Foam::cornerDetectionModels::fluxBased::detectCorners()
{
    if (!init_ || mesh().moving())
    {
        // Identify and store corner edges based on face normals
        cornerEdges_ = identifyCornerEdges();
        init_ = true;
    }


    // Identify and store corner faces based on edge flux sign
    this->setCornerFaces(identifyCornerFaces(cornerEdges_));


    // Calculate and store corner face angles
    const bitSet& cornerFaces = this->getCornerFaces();
    this->setCornerAngles
    (
        calcCornerAngles(cornerFaces, cornerEdges_)
    );


    // Write edges and faces as OBJ sets for debug purposes, if need be
    if (debug && mesh().time().writeTime())
    {
        writeEdgesAndFaces();
    }

    return true;
}


bool Foam::cornerDetectionModels::fluxBased::read(const dictionary& dict)
{
    if (!cornerDetectionModel::read(dict))
    {
        return false;
    }

    // Force the re-identification of corner edges/faces
    init_ = false;

    return true;
}


// ************************************************************************* //
