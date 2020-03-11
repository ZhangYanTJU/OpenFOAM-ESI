/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "averageNeighbourFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "cellAspectRatio.H"
#include "syncTools.H"
#include "polyMeshTools.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(averageNeighbourFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        averageNeighbourFvGeometryScheme,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::averageNeighbourFvGeometryScheme::makePyrHeights
(
    const pointField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceNormals,

    scalarField& ownHeight,
    scalarField& neiHeight
) const
{
    ownHeight.setSize(mesh_.nFaces());
    neiHeight.setSize(mesh_.nInternalFaces());

    typedef Vector<solveScalar> solveVector;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const solveVector n = faceNormals[facei];
        const solveVector fc = faceCentres[facei];
        ownHeight[facei] = ((fc-cellCentres[own[facei]])&n);
        neiHeight[facei] = ((cellCentres[nei[facei]]-fc)&n);
    }

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const solveVector n = faceNormals[facei];
        const solveVector fc = faceCentres[facei];
        ownHeight[facei] = ((fc-cellCentres[own[facei]])&n);
    }
}


Foam::label Foam::averageNeighbourFvGeometryScheme::clipPyramids
(
    const pointField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceNormals,

    const scalarField& minOwnHeight,
    const scalarField& minNeiHeight,

    vectorField& correction
) const
{
    // Clip correction vector if any pyramid becomes too small. Return number of
    // cells clipped

    typedef Vector<solveScalar> solveVector;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    label nClipped = 0;
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const vector& n = faceNormals[facei];
        const point& fc = faceCentres[facei];

        const label ownCelli = own[facei];
        if (correction[ownCelli] != vector::zero)
        {
            const solveVector ownCc(cellCentres[ownCelli]+correction[ownCelli]);
            const scalar ownHeight = ((fc-ownCc)&n);
            if (ownHeight < minOwnHeight[facei])
            {
                //Pout<< "    internalface:" << fc
                //    << " own:" << ownCc
                //    << " pyrHeight:" << ownHeight
                //    << " minHeight:" << minOwnHeight[facei]
                //    << endl;
                correction[ownCelli] = vector::zero;
                nClipped++;
            }
        }

        const label neiCelli = nei[facei];
        if (correction[neiCelli] != vector::zero)
        {
            const solveVector neiCc(cellCentres[neiCelli]+correction[neiCelli]);
            const scalar neiHeight = ((neiCc-fc)&n);
            if (neiHeight < minNeiHeight[facei])
            {
                //Pout<< "    internalface:" << fc
                //    << " nei:" << neiCc
                //    << " pyrHeight:" << neiHeight
                //    << " minHeight:" << minNeiHeight[facei]
                //    << endl;
                correction[neiCelli] = vector::zero;
                nClipped++;
            }
        }
    }

    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        const vector& n = faceNormals[facei];
        const point& fc = faceCentres[facei];

        const label ownCelli = own[facei];
        if (correction[ownCelli] != vector::zero)
        {
            const solveVector ownCc(cellCentres[ownCelli]+correction[ownCelli]);
            const scalar ownHeight = ((fc-ownCc)&n);
            if (ownHeight < minOwnHeight[facei])
            {
                //Pout<< "    boundaryface:" << fc
                //    << " own:" << ownCc
                //    << " pyrHeight:" << ownHeight
                //    << " minHeight:" << minOwnHeight[facei]
                //    << endl;
                correction[ownCelli] = vector::zero;
                nClipped++;
            }
        }
    }
    return returnReduce(nClipped, sumOp<label>());
}


Foam::tmp<Foam::pointField>
Foam::averageNeighbourFvGeometryScheme::averageNeighbourCentres
(
    const pointField& cellCentres,
    const vectorField& faceNormals,
    const scalarField& faceWeights
) const
{
    typedef Vector<solveScalar> solveVector;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    tmp<pointField> tcc(new pointField(mesh_.nCells(), Zero));
    pointField& cc = tcc.ref();

    Field<solveScalar> cellWeights(mesh_.nCells(), Zero);

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const vector& n = faceNormals[facei];
        const point& ownCc = cellCentres[own[facei]];
        const point& neiCc = cellCentres[nei[facei]];

        solveVector d(neiCc-ownCc);

        // 1. Normalise contribution. This increases actual non-ortho
        // since it does not 'see' the tangential offset of neighbours
        //neiCc = ownCc + (d&n)*n;

        // 2. Remove normal contribution, i.e. get tangential vector
        //    (= non-ortho correction vector?)
        d -= (d&n)*n;

        // Apply half to both sides (as a correction)
        // Note: should this be linear weights instead of 0.5?
        const scalar w = 0.5*faceWeights[facei];
        cc[own[facei]] += w*d;
        cellWeights[own[facei]] += w;

        cc[nei[facei]] -= w*d;
        cellWeights[nei[facei]] += w;
    }


    // Boundary faces. Bypass stored cell centres
    pointField neiCellCentres;
    syncTools::swapBoundaryCellPositions(mesh_, cellCentres, neiCellCentres);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            const labelUList& fc = pp.faceCells();

            forAll(fc, i)
            {
                const label meshFacei = pp.start()-mesh_.nInternalFaces()+i;
                const label bFacei = meshFacei-mesh_.nInternalFaces();

                const vector& n = faceNormals[meshFacei];

                const point& ownCc = cellCentres[fc[i]];
                const point& neiCc = neiCellCentres[bFacei];

                solveVector d(neiCc-ownCc);

                // 1. Normalise contribution. This increases actual non-ortho
                // since it does not 'see' the tangential offset of neighbours
                //neiCc = ownCc + (d&n)*n;

                // 2. Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                // Apply half to both sides (as a correction)
                const scalar w = 0.5*faceWeights[meshFacei];
                cc[fc[i]] += w*d;
                cellWeights[fc[i]] += w;
            }
        }
    }

    // Now cc is still the correction vector
    cc /= cellWeights;

    // Add correction vector to the input cell centres
    cc += cellCentres;

    return tcc;
}


Foam::tmp<Foam::pointField>
Foam::averageNeighbourFvGeometryScheme::averageCentres
(
    const pointField& cellCentres,
    const pointField& faceCentres,
    const vectorField& faceNormals
) const
{
    typedef Vector<solveScalar> solveVector;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    tmp<pointField> tnewFc(new pointField(faceCentres));
    pointField& newFc = tnewFc.ref();

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        const vector& n = faceNormals[facei];
        const point& oldFc = faceCentres[facei];

        const solveVector ownCc(cellCentres[own[facei]]);
        const solveVector neiCc(cellCentres[nei[facei]]);

        solveVector deltaCc(neiCc-ownCc);
        solveVector deltaFc(oldFc-ownCc);

        //solveVector d(neiCc-ownCc);
        //// 1. Normalise contribution. This increases actual non-ortho
        //// since it does not 'see' the tangential offset of neighbours
        ////neiCc = ownCc + s*n;
        //
        //// 2. Remove normal contribution, i.e. get tangential vector
        ////    (= non-ortho correction vector?)
        //d -= s*n;
        //newFc[facei] = faceCentres[facei]+d;

        // Get linear weight (normal distance to face)
        const solveScalar f = (deltaFc&n)/(deltaCc&n);
        const solveVector avgCc((1.0-f)*ownCc + f*neiCc);

        solveVector d(avgCc-oldFc);
        // Remove normal contribution, i.e. get tangential vector
        //    (= non-ortho correction vector?)
        d -= (d&n)*n;

        newFc[facei] = oldFc + d;
    }


    // Boundary faces. Bypass stored cell centres
    pointField neiCellCentres;
    syncTools::swapBoundaryCellPositions(mesh_, cellCentres, neiCellCentres);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelUList& fc = pp.faceCells();

        if (pp.coupled())
        {
            forAll(fc, i)
            {
                // Same as internal faces
                const label facei = pp.start()+i;
                const label bFacei = facei-mesh_.nInternalFaces();

                const vector& n = faceNormals[facei];
                const point& oldFc = faceCentres[facei];

                const solveVector ownCc(cellCentres[fc[i]]);
                const solveVector neiCc(neiCellCentres[bFacei]);

                solveVector deltaCc(neiCc-ownCc);
                solveVector deltaFc(oldFc-ownCc);

                // Get linear weight (normal distance to face)
                const solveScalar f = (deltaFc&n)/(deltaCc&n);
                const solveVector avgCc((1.0-f)*ownCc + f*neiCc);

                solveVector d(avgCc-oldFc);
                // Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                newFc[facei] = oldFc + d;
            }
        }
        else
        {
            // Zero-grad?
            forAll(fc, i)
            {
                const label facei = pp.start()+i;

                const vector& n = faceNormals[facei];
                const point& oldFc = faceCentres[facei];
                const solveVector ownCc(cellCentres[fc[i]]);

                solveVector d(ownCc-oldFc);
                // Remove normal contribution, i.e. get tangential vector
                //    (= non-ortho correction vector?)
                d -= (d&n)*n;

                newFc[facei] = oldFc+d;
            }
        }
    }

    return tnewFc;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::averageNeighbourFvGeometryScheme::averageNeighbourFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    highAspectRatioFvGeometryScheme(mesh, dict),
    nIters_(dict.getOrDefault<label>("nIters", 1)),
    relax_(dict.get<scalar>("relax"))
{
    // Force local calculation
    movePoints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::averageNeighbourFvGeometryScheme::movePoints()
{
    if (debug)
    {
        Pout<< "averageNeighbourFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    //if
    //(
    //   !mesh_.hasCellCentres()
    //&& !mesh_.hasFaceCentres()
    //&& !mesh_.hasCellVolumes()
    //&& !mesh_.hasFaceAreas()
    //)
    {
        highAspectRatioFvGeometryScheme::movePoints();

        // Note: at this point the highAspectRatioFvGeometryScheme constructor
        //       will have already reset the primitive geometry!

        vectorField faceAreas(mesh_.faceAreas());
        const scalarField magFaceAreas(mag(faceAreas));
        const vectorField faceNormals(faceAreas/magFaceAreas);


        // Calculate aspectratio weights
        // - 0 if aratio < minAspect_
        // - 1 if aratio >= maxAspect_
        scalarField cellWeight, faceWeight;
        calcAspectRatioWeights(cellWeight, faceWeight);

        // Relaxation
        cellWeight *= relax_;
        //faceWeight *= relax_;

        // Calculate current pyramid heights
        scalarField minOwnHeight;
        scalarField minNeiHeight;
        makePyrHeights
        (
            mesh_.cellCentres(),
            mesh_.faceCentres(),
            faceNormals,

            minOwnHeight,
            minNeiHeight
        );

        // How much is the cell centre to vary inside the cell.
        minOwnHeight *= 0.5;
        minNeiHeight *= 0.5;

        pointField cellCentres(mesh_.cellCentres());

        // Modify cell centres to be more in-line with the face normals
        pointField newCellCentres(mesh_.cellCentres());
        for (label iter = 0; iter < nIters_; iter++)
        {
            // Get neighbour average. Clip to limit change in pyramid height
            // (minOwnWeight, minNeiWeight)
            tmp<pointField> tcc
            (
                averageNeighbourCentres
                (
                    cellCentres,
                    faceNormals,
                    magFaceAreas
                )
            );

            // Calculate correction for cell centres. Leave low-aspect
            // ratio cells unaffected (i.e. correction = 0)
            vectorField correction(cellWeight*(tcc-cellCentres));

            // Clip correction vector if pyramid becomes too small
            const label nClipped = clipPyramids
            (
                cellCentres,
                mesh_.faceCentres(),
                faceNormals,

                minOwnHeight, // minimum owner pyramid height. Usually fraction
                minNeiHeight, // of starting mesh

                correction
            );
            //DebugVar(nClipped);

            cellCentres += correction;

            if (debug)
            {
                const scalarField magCorrection(mag(correction));
                const scalarField faceOrthogonality
                (
                    polyMeshTools::faceOrthogonality
                    (
                        mesh_,
                        faceAreas,
                        cellCentres
                    )
                );
                const scalarField nonOrthoAngle
                (
                    radToDeg
                    (
                        Foam::acos(min(scalar(1), faceOrthogonality))
                    )
                );
                Pout<< "    iter:" << iter
                    << " nClipped:" << nClipped
                    << " average displacement:" << gAverage(magCorrection)
                    << " non-ortho angle : average:" << gAverage(nonOrthoAngle)
                    << " max:" << gMax(nonOrthoAngle) << endl;
            }
        }

        tmp<pointField> tfc
        (
            averageCentres
            (
                cellCentres,
                mesh_.faceCentres(),
                faceNormals
            )
        );
        vectorField faceCentres
        (
            (1.0-faceWeight)*mesh_.faceCentres()
          + faceWeight*tfc
        );

        if (debug)
        {
            Pout<< "averageNeighbourFvGeometryScheme::movePoints() :"
                << " averageNeighbour weight"
                << " max:" << gMax(cellWeight) << " min:" << gMin(cellWeight)
                << " average:" << gAverage(cellWeight) << endl;

            // Dump lines from old to new location
            const fileName tp(mesh_.time().timePath());
            mkDir(tp);
            OBJstream str(tp/"averageNeighbourCellCentres.obj");
            Pout<< "Writing lines from old to new cell centre to " << str.name()
                << endl;
            forAll(mesh_.cellCentres(), celli)
            {
                const point& oldCc = mesh_.cellCentres()[celli];
                const point& newCc = cellCentres[celli];
                str.write(linePointRef(oldCc, newCc));
            }
        }
        if (debug)
        {
            // Dump lines from old to new location
            const fileName tp(mesh_.time().timePath());
            OBJstream str(tp/"averageFaceCentres.obj");
            Pout<< "Writing lines from old to new face centre to " << str.name()
                << endl;
            forAll(mesh_.faceCentres(), facei)
            {
                const point& oldFc = mesh_.faceCentres()[facei];
                const point& newFc = faceCentres[facei];
                str.write(linePointRef(oldFc, newFc));
            }
        }

        scalarField cellVolumes(mesh_.cellVolumes());

        // Store on primitiveMesh
        //const_cast<fvMesh&>(mesh_).clearGeom();
        const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
        (
            std::move(faceCentres),
            std::move(faceAreas),
            std::move(cellCentres),
            std::move(cellVolumes)
        );
    }
}


// ************************************************************************* //
