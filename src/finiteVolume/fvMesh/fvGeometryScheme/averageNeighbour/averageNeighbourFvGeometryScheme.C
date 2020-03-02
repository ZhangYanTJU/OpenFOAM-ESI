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

Foam::tmp<Foam::pointField>
Foam::averageNeighbourFvGeometryScheme::averageNeighbourCentres
(
    const pointField& cellCentres,
    const vectorField& faceNormals,
    const scalarField& faceWeights
) const
{
    if (debug)
    {
        Pout<< "highAspectRatioFvGeometryScheme::averageNeighbourCentres() : "
            << "calculating weighted neighbouring cell centre" << endl;
    }

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

        const solveVector myCc(cellCentres[own[facei]]);
        const solveVector nbrCc(cellCentres[nei[facei]]);

        solveVector d(nbrCc-myCc);

        // 1. Normalise contribution. This increases actual non-ortho
        // since it does not 'see' the tangential offset of neighbours
        //nbrCc = myCc + (d&n)*n;

        // 2. Remove normal contribution, i.e. get tangential vector
        //    (= non-ortho correction vector?)
        d -= (d&n)*n;

        // Apply half to both sides (as a correction)
        const scalar w = 0.5*faceWeights[facei];
        cc[own[facei]] += w*d;
        cellWeights[own[facei]] += w;

        cc[nei[facei]] -= w*d;
        cellWeights[nei[facei]] += w;
    }


    // Boundary faces. Bypass stored cell centres
    pointField nbrCellCentres;
    syncTools::swapBoundaryCellPositions(mesh_, cellCentres, nbrCellCentres);

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

                const solveVector myCc(cellCentres[fc[i]]);
                const solveVector nbrCc(nbrCellCentres[bFacei]);

                solveVector d(nbrCc-myCc);

                // 1. Normalise contribution. This increases actual non-ortho
                // since it does not 'see' the tangential offset of neighbours
                //nbrCc = myCc + (d&n)*n;

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

    cc /= cellWeights;
    cc += cellCentres;
    return tcc;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::averageNeighbourFvGeometryScheme::averageNeighbourFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    highAspectRatioFvGeometryScheme(mesh, dict)
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

    if
    (
       !mesh_.hasCellCentres()
    && !mesh_.hasFaceCentres()
    && !mesh_.hasCellVolumes()
    && !mesh_.hasFaceAreas()
    )
    {
        vectorField faceAreas(mesh_.faceAreas());
        scalarField cellVolumes(mesh_.cellVolumes());
        const scalarField magFaceAreas(mag(faceAreas));

        // Calculate weighted average centroids
        pointField avgFaceCentres;
        pointField avgCellCentres;
        makeAverageCentres
        (
            mesh_,
            mesh_.points(),
            faceAreas,
            magFaceAreas,
            avgFaceCentres,
            avgCellCentres
        );

        // Modify cell centres to be more in-line with the face normals
        tmp<pointField> tcc
        (
            averageNeighbourCentres
            (
                avgCellCentres,
                faceAreas/magFaceAreas,
                magFaceAreas
            )
        );


        const cellAspectRatio aratio(mesh_);

        // Weighting for correction
        // - 0 if aratio < minAspect_
        // - 1 if aratio >= maxAspect_
        const scalarField cellWeight
        (
            max
            (
                min
                (
                    (aratio-minAspect_)/(maxAspect_-minAspect_),
                    1.0
                ),
                0.0
            )
        );

        scalarField faceWeight(mesh_.nFaces());
        {
            for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
            {
                const label own = mesh_.faceOwner()[facei];
                const label nei = mesh_.faceNeighbour()[facei];
                faceWeight[facei] = max(cellWeight[own], cellWeight[nei]);
            }
            scalarField nbrCellWeight;
            syncTools::swapBoundaryCellList
            (
                mesh_,
                cellWeight,
                nbrCellWeight
            );
            for
            (
                label facei = mesh_.nInternalFaces();
                facei < mesh_.nFaces();
                facei++
            )
            {
                const label own = mesh_.faceOwner()[facei];
                const label bFacei = facei-mesh_.nInternalFaces();
                faceWeight[facei] = max(cellWeight[own], nbrCellWeight[bFacei]);
            }
        }


        // Weight with average ones

        vectorField faceCentres
        (
            (1.0-faceWeight)*mesh_.faceCentres()
          + faceWeight*avgFaceCentres
        );
        vectorField cellCentres
        (
            (1.0-cellWeight)*mesh_.cellCentres()
          + cellWeight*tcc
        );


        if (debug)
        {
            Pout<< "averageNeighbourFvGeometryScheme::movePoints() :"
                << " averageNeighbour weight"
                << " max:" << gMax(cellWeight) << " min:" << gMin(cellWeight)
                << " average:" << gAverage(cellWeight) << endl;
        }

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
