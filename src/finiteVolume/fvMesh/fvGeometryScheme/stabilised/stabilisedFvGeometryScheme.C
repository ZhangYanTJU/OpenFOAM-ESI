/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
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

#include "stabilisedFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "PrecisionAdaptor.H"
#include "primitiveMeshTools.H"
#include "triangle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stabilisedFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        stabilisedFvGeometryScheme,
        dict
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Average of points
// Note: use double precision to avoid overflows when summing
static inline Foam::doubleVector pointsAverage
(
    const UList<point>& points,
    const labelUList& pointLabels
)
{
    doubleVector avg(Zero);

    if (const auto n = pointLabels.size(); n)
    {
        for (const auto pointi : pointLabels)
        {
            avg += points[pointi];
        }
        avg /= n;
    }

    return avg;
}

} // End namespace Foam



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::stabilisedFvGeometryScheme::makeFaceCentresAndAreas
(
    const polyMesh& mesh,
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas
)
{
    const faceList& fs = mesh.faces();

    forAll(fs, facei)
    {
        const auto& f = fs[facei];
        const label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            const triPointRef tri(p, f[0], f[1], f[2]);
            fCtrs[facei] = tri.centre();
            fAreas[facei] = tri.areaNormal();
        }

        // For more complex faces, decompose into triangles
        else
        {
            // Estimated centre by averaging the face points
            const solveVector fCentre(pointsAverage(p, f));

            // Compute the face area normal and unit normal by summing up the
            // normals of the triangles formed by connecting each edge to the
            // point average.
            solveVector sumA(Zero);
            for (label pi = 0; pi < nPoints; pi++)
            {
                const label nextPi(pi == nPoints-1 ? 0 : pi+1);
                const solveVector nextPoint(p[f[nextPi]]);
                const solveVector thisPoint(p[f[pi]]);

                const solveVector a =
                    (nextPoint - thisPoint)^(fCentre - thisPoint);

                sumA += a;
            }
            const solveVector sumAHat = normalised(sumA);

            // Compute the area-weighted sum of the triangle centres. Note use
            // the triangle area projected in the direction of the face normal
            // as the weight, *not* the triangle area magnitude. Only the
            // former makes the calculation independent of the initial estimate.
            solveScalar sumAn(0);
            solveVector sumAnc(Zero);
            for (label pi = 0; pi < nPoints; pi++)
            {
                const label nextPi(pi == nPoints-1 ? 0 : pi+1);
                const solveVector nextPoint(p[f[nextPi]]);
                const solveVector thisPoint(p[f[pi]]);

                const solveVector c = thisPoint + nextPoint + fCentre;
                const solveVector a =
                    (nextPoint - thisPoint)^(fCentre - thisPoint);

                const scalar an = a & sumAHat;

                sumAn += an;
                sumAnc += an*c;
            }

            // Complete calculating centres and areas. If the face is too small
            // for the sums to be reliably divided then just set the centre to
            // the initial estimate.
            if (sumAn > ROOTVSMALL)
            {
                fCtrs[facei] = (1.0/3.0)*sumAnc/sumAn;
            }
            else
            {
                fCtrs[facei] = fCentre;
            }
            fAreas[facei] = 0.5*sumA;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stabilisedFvGeometryScheme::stabilisedFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicFvGeometryScheme(mesh, dict)
{
    // Force local calculation
    movePoints();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stabilisedFvGeometryScheme::movePoints()
{
    fvGeometryScheme::movePoints();

    if (debug)
    {
        Pout<< "stabilisedFvGeometryScheme::movePoints() : "
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
        vectorField faceCentres(mesh_.nFaces());
        vectorField faceAreas(mesh_.nFaces());

        makeFaceCentresAndAreas
        (
            mesh_,
            mesh_.points(),
            faceCentres,
            faceAreas
        );

        vectorField cellCentres(mesh_.nCells());
        scalarField cellVolumes(mesh_.nCells());

        primitiveMeshTools::makeCellCentresAndVols
        (
            mesh_,
            faceCentres,
            faceAreas,
            cellCentres,
            cellVolumes
        );

        const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
        (
            std::move(faceCentres),
            std::move(faceAreas),
            std::move(cellCentres),
            std::move(cellVolumes)
        );
    }
}


bool Foam::stabilisedFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    makeFaceCentresAndAreas
    (
        mesh_,
        points,
        faceCentres,
        faceAreas
    );

    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh_,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes
    );

    return true;
}


// ************************************************************************* //
