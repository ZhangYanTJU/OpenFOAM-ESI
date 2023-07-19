/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "topOVariablesBase.H"
#include "pointMesh.H"
#include "topOZones.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "coupledFvPatch.H"
#include "emptyFvPatch.H"
#include "MeshedSurfaceProxy.H"
#include "surfaceWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topOVariablesBase, 0);
}


// * * * * * * * * * * * Protected Members Functions * * * * * * * * * * * * //

bool Foam::topOVariablesBase::addCutBoundaryFaceToIsoline
(
    const label facei,
    const cutFaceIso& cutFace,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs,
    label& nIsoSurfPts
) const
{
    const fvMesh& mesh = zones_.mesh();
    // If this is a boundary face being cut, append it to the iso-surface
    if (!mesh.isInternalFace(facei))
    {
        const label cutPatchi = mesh.boundaryMesh().whichPatch(facei);
        const fvPatch& cutPatch = mesh.boundary()[cutPatchi];
        if
        (
            !isA<coupledFvPatch>(cutPatch)
         && !isA<emptyFvPatch>(cutPatch)
         && (cutFace.subFacePoints().size() > 2)
        )
        {
            const DynamicList<point>& facePts = cutFace.subFacePoints();
            face isoFace(identity(facePts.size(), nIsoSurfPts));
            isoSurfPts.append(facePts);
            isoSurfFaces.append(isoFace);
            zoneIDs.append(cutPatchi);
            nIsoSurfPts += facePts.size();
            return true;
        }
    }
    return false;
}


bool Foam::topOVariablesBase::addCuttingFaceToIsoline
(
    const DynamicList<point>& facePoints,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs,
    label& nIsoSurfPts,
    const label nSerialPatches
) const
{
    if (facePoints.size() > 2)
    {
        face isoFace(identity(facePoints.size(), nIsoSurfPts));
        isoSurfPts.append(facePoints);
        isoSurfFaces.append(isoFace);
        zoneIDs.append(nSerialPatches);
        nIsoSurfPts += facePoints.size();
        return true;
    }
    return false;
}


void Foam::topOVariablesBase::addBoundaryFacesToIsoline
(
    const pointScalarField& pointY,
    const Map<label>& addedFaces,
    const scalar isoValue,
    DynamicList<vector>& isoSurfPts,
    DynamicList<face>& isoSurfFaces,
    DynamicList<label>& zoneIDs,
    label& nIsoSurfPts
) const
{
    const fvMesh& mesh = zones_.mesh();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if
        (
            !isA<emptyFvPatch>(patch)
         && !isA<coupledFvPatch>(patch)
        )
        {
            const label start = patch.start();
            forAll(patch, facei)
            {
                const label gFacei = start + facei;
                // On rare occasions, a face value might be negative but all
                // point values are positive. In such a case, the face won't be
                // cut and is seen as a fluid face. Hence, the point values are
                // used to determine the status of the boundary face.
                bool isFluid(true);
                for (const label pti : mesh.faces()[gFacei])
                {
                    isFluid = isFluid && (pointY[pti] >= isoValue);
                }
                if (isFluid && !addedFaces.found(gFacei))
                {
                    const pointField facePoints = faces[gFacei].points(points);
                    const face isoFace
                        (identity(facePoints.size(), nIsoSurfPts));
                    isoSurfPts.append(facePoints);
                    isoSurfFaces.append(isoFace);
                    zoneIDs.append(patchi);
                    nIsoSurfPts += facePoints.size();
                }
            }
        }
    }
}


void Foam::topOVariablesBase::writeSurfaceFiles
(
    DynamicList<vector>& patchPoints,
    DynamicList<face>& patchFaces,
    DynamicList<label>& zoneIDs,
    const label nSerialPatches
) const
{
    const fvMesh& mesh = zones_.mesh();

    pointField pts(std::move(patchPoints));
    faceList faces(std::move(patchFaces));
    labelList zoneIds(std::move(zoneIDs));

    // Write vtp file
    const word timeName = mesh.time().timeName();
    fileName localName("topOIsoSurface" + timeName);
    fileName fname(isoSurfFolder_/localName);
    vtk::surfaceWriter vtkFile(pts, faces, fname);
    vtkFile.beginFile();
    vtkFile.writeGeometry();
    vtkFile.beginCellData(1);
    vtkFile.writeCellData("zoneIds", zoneIds);
    vtkFile.close();

    // It appears that there is no interface to write a multi-zone stl file (?)
    // Follow a process similar to proxySurfaceWriter, but also use the zoneIds.

    // Dummy faceIds.
    // faceMap is built after merging the geometry from all processors, to
    // be based on the global addressing
    labelList faceIds;

    autoPtr<meshedSurf> surf(nullptr);
    if (Pstream::parRun())
    {
        surf.reset
        (
            autoPtr<meshedSurf>::NewFrom<mergedSurf>
            (
                pts, faces, zoneIds, faceIds, surfaceWriter::defaultMergeDim
            )
       );
    }
    else
    {
        surf.reset
        (
            autoPtr<meshedSurf>::NewFrom<meshedSurfRef>
            (
                pts, faces, zoneIds, faceIds
            )
        );
    }

    if (Pstream::master())
    {
        const faceList& surfFaces = surf().faces();
        const labelList& surfZoneIds = surf().zoneIds();

        // Size per zone
        labelList zoneSizes(nSerialPatches + 1, 0);
        forAll(surfFaces, fi)
        {
            ++zoneSizes[surfZoneIds[fi]];
        }

        // Faces passed in previous zones
        labelList cumulZoneSizes(surfZoneIds.size(), 0);
        for (label zi = 1; zi < surfZoneIds.size() - 1; ++zi)
        {
            cumulZoneSizes[zi] = cumulZoneSizes[zi - 1] + zoneSizes[zi - 1];
        }

        // Construction of the faceMap
        // ---------------------------
        labelList faceMap(surfFaces.size(), -1);
        // Faces visited so far per zone
        labelList passedFacesPerZone(surfZoneIds.size(), 0);

        forAll(surfFaces, fi)
        {
            const label zi = surfZoneIds[fi];
            faceMap[cumulZoneSizes[zi] + passedFacesPerZone[zi]] = fi;
            ++passedFacesPerZone[zi];
        }

        // Construction of the surfZoneList
        // --------------------------------
        List<surfZone> surfZones(nSerialPatches + 1);
        label zi = 0;
        forAll(mesh.boundary(), pi)
        {
            const fvPatch& patch = mesh.boundary()[pi];
            const word& name = patch.name();
            if (!isA<coupledFvPatch>(patch) && !isA<emptyFvPatch>(patch))
            {
                surfZones[zi] =
                    surfZone(name, zoneSizes[zi], cumulZoneSizes[zi], zi);
                zi++;
            }
        }
        surfZones[zi] =
            surfZone
            (
                "topOPatch",
                zoneSizes[nSerialPatches],
                cumulZoneSizes[nSerialPatches],
                zi
            );
        surfZones.setSize(zi + 1);

        MeshedSurfaceProxy<face>
        (
            surf().points(),
            surfFaces,
            surfZones,
            faceMap
        ).write
        (
            fname + ".stl"
            // BINARY seems to loose the patch names
            //IOstreamOption
            //(
            //    IOstreamOption::streamFormat::ASCII,
            //    IOstreamOption::compressionType::COMPRESSED
            //)
        );
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topOVariablesBase::topOVariablesBase
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    localIOdictionary
    (
        IOobject
        (
            "topoVars",
            mesh.time().timeName(),
            fileName("uniform"),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        word::null
    ),
    zones_(mesh, dict),
    isoSurfFolder_
        (mesh.time().globalPath()/"optimisation"/"topOIsoSurfaces")
{
    mkDir(isoSurfFolder_);
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::topOVariablesBase::sourceTerm
(
    DimensionedField<scalar, volMesh>& field,
    const topOInterpolationFunction& interpolationFunc,
    const scalar betaMax,
    const word& interpolationFieldName
) const
{
    // Interpolate based on beta values
    const scalarField& beta = this->beta().primitiveField();
    interpolationFunc.interpolate(beta, field.field());
    field *= betaMax;
}


void Foam::topOVariablesBase::sourceTermSensitivities
(
    scalarField& sens,
    const topOInterpolationFunction& interpolationFunc,
    const scalar betaMax,
    const word& designVariablesName,
    const word& interpolationFieldName
) const
{
    if (designVariablesName == "beta")
    {
        // Multiply with derivative of the interpolation function
        const scalarField& beta = this->beta().primitiveField();
        sens *= betaMax*interpolationFunc.derivative(beta);
    }
}


void Foam::topOVariablesBase::writeFluidSolidInterface
(
    const volScalarField& indicator,
    const scalar isoValue,
    labelList& changedFaces,
    List<wallPoint>& changedFacesInfo
)
{
    const fvMesh& mesh = zones_.mesh();

    // Current number of interface faces
    label nSurfFaces(0);

    // Indicator at the mesh points, used to compute the iso-surface
    pointScalarField pointY(volPointInterpolation(mesh).interpolate(indicator));
    if (debug && mesh.time().writeTime())
    {
        pointY.rename("pointY");
        pointY.write();
    }

    // Storage for the cutting face/cell mechanism
    cutCellIso cutCell(mesh, pointY);
    cutFaceIso cutFace(mesh, pointY);
    DynamicList<face> isoSurfFaces(mesh.nFaces()/100);
    DynamicList<point> isoSurfPts(mesh.nPoints()/100);
    DynamicList<label> zoneIDs(mesh.nFaces()/100);
    label nIsoSurfPts(0);

    // Number of patches before processor ones
    label nSerialPatches = mesh.boundaryMesh().nNonProcessor();

    // Map between iso-surface and mesh faces (internal and boundary)
    Map<label> addedFaces;

    // Loop over cells and check whether they are cut by zero iso-surface
    const cellList& cells = mesh.cells();
    forAll(cells, celli)
    {
        if (!cutCell.calcSubCell(celli, isoValue))
        {
            const vector& cuttingCf = cutCell.faceCentre();
            vector cuttingNf = cutCell.faceArea();
            cuttingNf.normalise();

            // Loop over faces of the cut cell and check whether they are cut
            // too, setting their distance from the cutting face
            for (const label facei : cells[celli])
            {
                if (!cutFace.calcSubFace(facei, isoValue))
                {
                    const vector& Cf = mesh.Cf()[facei];
                    const scalar distSqr =
                        magSqr((Cf - cuttingCf) & cuttingNf);
                    const vector lsPoint =
                        Cf + ((cuttingCf - Cf) & cuttingNf)*cuttingNf;

                    if (!addedFaces.found(facei))
                    {
                        // If the face being cut is a boundary one, part of
                        // it belongs to the iso-surface
                        bool addedToIsoSurf = addCutBoundaryFaceToIsoline
                        (
                            facei,
                            cutFace,
                            isoSurfPts,
                            isoSurfFaces,
                            zoneIDs,
                            nIsoSurfPts
                        );
                        if (mesh.isInternalFace(facei) || addedToIsoSurf)
                        {
                            addedFaces.insert(facei, nSurfFaces);
                        }
                        changedFacesInfo[nSurfFaces] =
                            wallPoint(lsPoint, distSqr);
                        changedFaces[nSurfFaces++] = facei;
                    }
                    else
                    {
                        const label visitedFace = addedFaces.at(facei);
                        if (distSqr < changedFacesInfo[visitedFace].distSqr())
                        {
                            changedFacesInfo[visitedFace] =
                                wallPoint(lsPoint, distSqr);
                        }
                    }
                }
            }

            addCuttingFaceToIsoline
            (
                cutCell.facePoints(),
                isoSurfPts,
                isoSurfFaces,
                zoneIDs,
                nIsoSurfPts,
                nSerialPatches
            );

        }
    }

    // Insert wall faces if they belong to the outer boundary
    labelHashSet wallPatchIDs =
        mesh.boundaryMesh().findPatchIDs<wallPolyPatch>();
    for (const label patchi : wallPatchIDs)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const labelList& faceCells = patch.faceCells();
        const label start = patch.start();
        const vectorField& Cf = patch.Cf();

        forAll(Cf, facei)
        {
            if
            (
                indicator[faceCells[facei]] >= 0 &&
               !addedFaces.found(start + facei)
            )
            {
                changedFaces[nSurfFaces] = start + facei;
                changedFacesInfo[nSurfFaces] = wallPoint(Cf[facei], 0);
                ++nSurfFaces;
            }
        }
    }

    // Add boundary faces on the initial domain with a positive values on all
    // points to the zero iso-surface
    addBoundaryFacesToIsoline
    (
        pointY,
        addedFaces,
        isoValue,
        isoSurfPts,
        isoSurfFaces,
        zoneIDs,
        nIsoSurfPts
    );

    changedFaces.setSize(nSurfFaces);
    changedFacesInfo.setSize(nSurfFaces);

    // vtp and multi-region stl files holding the current geometry
    writeSurfaceFiles(isoSurfPts, isoSurfFaces, zoneIDs, nSerialPatches);
}


// ************************************************************************* //
