/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

Application
    surfaceMeshExtract

Group
    grpSurfaceUtilities

Description
    Extract patch or faceZone surfaces from a polyMesh.
    Depending on output surface format triangulates faces.

    Region numbers on faces not guaranteed to be the same as the patch indices.

    Optionally only extracts named patches.

    Optionally filters out points on feature-edges and generates pointPatches
    for these - written to constant/pointMesh/boundary.

    If run in parallel, processor patches get filtered out by default and
    the mesh is merged (based on topology).

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "ListListOps.H"
#include "indirectPrimitivePatch.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "timeSelector.H"
#include "meshPointPatch.H"
#include "unitConversion.H"
#include "dummyTransform.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const wordRes& allow,
    const wordRes& deny
)
{
    // Name-based selection
    labelList indices
    (
        stringListOps::findMatching
        (
            patches,
            allow,
            deny,
            nameOp<polyPatch>()
        )
    );


    // Remove undesirable patches

    label count = 0;
    for (const label patchi : indices)
    {
        const polyPatch& pp = patches[patchi];

        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (Pstream::parRun() && bool(isA<processorPolyPatch>(pp)))
        {
            break; // No processor patches for parallel output
        }

        indices[count] = patchi;
        ++count;
    }

    indices.resize(count);

    return indices;
}


label addMeshPointPatches
(
    const polyMesh& mesh,
    const pointMesh& pMesh,

    const uindirectPrimitivePatch& allBoundary,
    const labelUList& faceToZone,
    const surfZoneIdentifierList& surfZones,

    const scalar edgeFeatureAngle,
    const scalar pointFeatureAngle,
    const bool verbose = true
)
{
    const pointBoundaryMesh& pointBm = pMesh.boundary();
    const label nPointPatches = pointBm.size();

    const globalMeshData& globalData = mesh.globalData();
    const indirectPrimitivePatch& cpp = globalData.coupledPatch();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const auto& mp = allBoundary.meshPoints();


    const vector nullVector(vector::uniform(0));
    const auto assignNonNull = [&](vector& x, const vector& y)
    {
        if (x == nullVector && y != nullVector)
        {
            x = y;
        }
    };

    // Calculate parallel-consistent point normals (as unweighted average
    // of faceNormals). Note: only valid on patch points, not on mesh points
    // that are coupled to these.
    const pointField pointNormals
    (
        PatchTools::pointNormals
        (
            mesh,
            allBoundary
        )
    );
    // Expand to all coupled points
    pointField meshPointNormals(mesh.nPoints(), nullVector);
    UIndirectList<vector>(meshPointNormals, mp) = pointNormals;
    syncTools::syncPointList
    (
        mesh,
        meshPointNormals,
        assignNonNull,
        nullVector
    );


    // Find correspondence between allBoundary and coupled edges
    labelList allEdges;
    labelList coupledEdges;
    bitSet sameEdgeOrientation;
    PatchTools::matchEdges
    (
        allBoundary,
        cpp,
        allEdges,
        coupledEdges,
        sameEdgeOrientation
    );



    // To construct the patches we need to know per edge
    //  - patch on either side (if topological feature edge)
    //  - faceNormal on either side (if feature angle)
    // We need to know per point
    //  - patches on all connected faces
    //  - faceNormals on all connected faces? And compare to average?
    //    or edge normals on all connected edges


    typedef Tuple2<label, vector> PN;
    const PN nullPN(-1, nullVector);



    // Point-based analysis
    // ~~~~~~~~~~~~~~~~~~~~


    // Collect per (mesh)point the zones (1, 2 or >2). Note: per mesh to
    // make it easier to sync. See edge-based code below where we explicitly
    // have to transfer from patch-edge to mesh-point etc. Note sure which one
    // fits better....
    labelPairList pointToZones(mesh.nPoints(), labelPair(-1, -1));
    {
        // Combine zones.
        const auto combineZones = [&](labelPair& x, const labelPair& y)
        {
            if (x == labelPair(-2, -2))
            {
                // Already marked
            }
            else if (y == labelPair(-2, -2))
            {
                x = y;
            }
            else
            {
                // Find first free slot
                if (x[0] == -1)
                {
                    if (y[0] != -1)
                    {
                        x[0] = y[0];
                    }
                    else
                    {
                        x[0] = y[1];
                    }
                }
                else if (x[1] == -1)
                {
                    if (y[0] != -1 && y[0] != x[0])
                    {
                        x[1] = y[0];
                    }
                    else if (y[1] != -1 && y[1] != x[1])
                    {
                        x[1] = y[1];
                    }
                }
                else
                {
                    // Both x slots filled. See if y adds a 3rd element
                    if (y[0] != -1 && y[0] != x[0] && y[0] != x[1])
                    {
                        x = labelPair(-2, -2);
                    }
                    else if (y[1] != -1 && y[1] != x[0] && y[1] != x[1])
                    {
                        x = labelPair(-2, -2);
                    }
                }
            }
        };


        forAll(allBoundary, facei)
        {
            const auto& f = allBoundary[facei];
            const label zonei = faceToZone[facei];
            for (const label pointi : f)
            {
                auto& pZones = pointToZones[pointi];

                if (pZones != labelPair(-2, -2) && !pZones.contains(zonei))
                {
                    if (pZones.first() == -1)
                    {
                        pZones.first() = zonei;
                    }
                    else if (pZones.second() == -1)
                    {
                        pZones.second() = zonei;
                    }
                    else
                    {
                        // Mark as >2 zones
                        pZones = labelPair(-2, -2);
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointToZones,
            combineZones,
            labelPair(-1, -1),
            dummyTransform()
        );
    }




    // Edge-based analysis
    // ~~~~~~~~~~~~~~~~~~~~

    // 1. Local analysis

    List<Pair<PN>> allEdgeToFaces
    (
        allBoundary.nEdges(),
        Pair<PN>(nullPN, nullPN)
    );
    {
        const auto& edgeFaces = allBoundary.edgeFaces();
        const auto& faceNormals = allBoundary.faceNormals();

        forAll(edgeFaces, edgei)
        {
            const auto& eFaces = edgeFaces[edgei];
            const vector& n0 = faceNormals[eFaces[0]];
            const label zone0 = faceToZone[eFaces[0]];
            if (eFaces.size() == 1)
            {
                allEdgeToFaces[edgei] = Pair<PN>(PN(zone0, n0), nullPN);
            }
            else
            {
                const vector& n1 = faceNormals[eFaces[1]];
                const label zone1 = faceToZone[eFaces[1]];
                allEdgeToFaces[edgei] = Pair<PN>
                (
                    PN(zone0, n0),
                    PN(zone1, n1)
                );
            }
        }
    }


    // 2. Sync across coupled patches

    {
        // Combine pair of normals
        const auto vectorPairMax = [&](Pair<PN>& x, const Pair<PN>& y)
        {
            if (x[0] == nullPN)
            {
                if (y[0] != nullPN)
                {
                    x[0] = y[0];
                }
                else
                {
                    x[0] = y[1];
                }
            }
            else if (x[1] == nullPN)
            {
                if (y[0] != nullPN && y[0] != x[0])
                {
                    x[1] = y[0];
                }
                else
                {
                    x[1] = y[1];
                }
            }
        };

        List<Pair<PN>> cppEdgeData
        (
            map.constructSize(),
            Pair<PN>(nullPN, nullPN)
        );
        UIndirectList<Pair<PN>>(cppEdgeData, coupledEdges) =
            UIndirectList<Pair<PN>>(allEdgeToFaces, allEdges);

        globalData.syncData
        (
            cppEdgeData,
            globalData.globalEdgeSlaves(),
            globalData.globalEdgeTransformedSlaves(),
            map,
            globalData.globalTransforms(),
            vectorPairMax,
            dummyTransform()
        );

        UIndirectList<Pair<PN>>(allEdgeToFaces, allEdges) =
            UIndirectList<Pair<PN>>(cppEdgeData, coupledEdges);
    }


    // Now we have all the per-patch edge information
    // - do inter-patch edges
    // - do feature-angle edges
    // Store on mesh points

    const auto assignNonNullPN = [&](PN& x, const PN& y)
    {
        if (x.second() == nullVector && y.second() != nullVector)
        {
            x = y;
        }
    };


    // Storing the normal for points that are on inter-patch edges
    vectorField patchEdgeNormal(mesh.nPoints(), nullVector);
    // Storing the normal for points that are on patch-internal feat edges
    List<PN> featEdgeNormal(mesh.nPoints(), nullPN);
    forAll(allEdgeToFaces, edgei)
    {
        const edge& e = allBoundary.edges()[edgei];
        const label mp0 = mp[e[0]];
        const label mp1 = mp[e[1]];

        const Pair<PN>& facesInfo = allEdgeToFaces[edgei];

        if (facesInfo[1] == nullPN)
        {
            // Boundary edge
            patchEdgeNormal[mp0] = pointNormals[e[0]];
            patchEdgeNormal[mp1] = pointNormals[e[1]];
        }
        else
        {
            if (facesInfo[0].first() != facesInfo[1].first())
            {
                // Inter-patch
                patchEdgeNormal[mp0] = pointNormals[e[0]];
                patchEdgeNormal[mp1] = pointNormals[e[1]];
            }
            else
            {
                // Same patch - check for feature angle

                const vector& n0 = facesInfo[0].second();
                const vector& n1 = facesInfo[1].second();

                if ((n0 & n1) < Foam::cos(degToRad(edgeFeatureAngle)))
                {
                    if (patchEdgeNormal[mp0] == nullVector)
                    {
                        featEdgeNormal[mp0] = PN
                        (
                            facesInfo[0].first(),   // zone
                            pointNormals[e[0]]
                        );
                    }
                    if (patchEdgeNormal[mp1] == nullVector)
                    {
                        featEdgeNormal[mp1] = PN
                        (
                            facesInfo[0].first(),   // zone
                            pointNormals[e[1]]
                        );
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        patchEdgeNormal,
        assignNonNull,
        nullVector
    );
    syncTools::syncPointList
    (
        mesh,
        featEdgeNormal,
        assignNonNullPN,
        nullPN,
        dummyTransform()
    );

    // Make sure that inter-patch points are not also in feature-edge
    // points. Note: not absolutely nessecary since all inter-patch points
    // will also be in the 'normal' facePointPatches.

    DynamicList<label> multiZoneMeshPoints(allBoundary.nPoints());
    forAll(pointToZones, pointi)
    {
        if (pointToZones[pointi] == labelPair(-2, -2))
        {
            multiZoneMeshPoints.append(pointi);
            // Unmark as feature angle point
            patchEdgeNormal[pointi] = nullVector;
            featEdgeNormal[pointi] = nullPN;
        }
    }


    DynamicList<label> twoZoneMeshPoints(allBoundary.nPoints());
    forAll(patchEdgeNormal, pointi)
    {
        if (patchEdgeNormal[pointi] != nullVector)
        {
            twoZoneMeshPoints.append(pointi);
            // Unmark as feature angle point
            featEdgeNormal[pointi] = nullPN;
        }
    }


    // Sort featEdgeNormal according to zone
    List<List<label>> zoneToMeshPoints(surfZones.size());
    List<vectorField> zoneToNormal(surfZones.size());
    {
        labelList sizes(surfZones.size(), 0);
        forAll(featEdgeNormal, pointi)
        {
            const auto& pInfo = featEdgeNormal[pointi];
            if (pInfo != nullPN)
            {
                const label zonei = pInfo.first();
                sizes[zonei]++;
            }
        }
        forAll(zoneToMeshPoints, zonei)
        {
            zoneToMeshPoints[zonei].setSize(sizes[zonei]);
            zoneToNormal[zonei].setSize(sizes[zonei]);
        }
        sizes = 0;
        forAll(featEdgeNormal, pointi)
        {
            const auto& pInfo = featEdgeNormal[pointi];
            if (pInfo != nullPN)
            {
                const label zonei = pInfo.first();
                const label index = sizes[zonei]++;
                zoneToMeshPoints[zonei][index] = pointi;
                zoneToNormal[zonei][index] = pInfo.second();
            }
        }
    }

    // Add patches
    forAll(zoneToMeshPoints, zonei)
    {
        const label nPoints =
            returnReduce(zoneToMeshPoints[zonei].size(), sumOp<label>());
        if (nPoints)
        {
            const word patchName(surfZones[zonei].name() + "Edges");

            const_cast<pointBoundaryMesh&>(pointBm).push_back
            (
                new meshPointPatch
                (
                    patchName,
                    zoneToMeshPoints[zonei],
                    zoneToNormal[zonei],
                    pointBm.size(),
                    pointBm,
                    meshPointPatch::typeName
                )
            );

            if (verbose)
            {
                Info<< "Added feature-edges pointPatch "
                    << pointBm.last().name()
                    << " index " << pointBm.last().index()
                    << " with " << nPoints << " points"
                    << endl;
            }
        }
    }


    // Add inter-patch points

    const word allPointPatchName("boundaryPoints");

    const label nMultiPoints =
        returnReduce(multiZoneMeshPoints.size(), sumOp<label>());
    if (nMultiPoints)
    {
        const_cast<pointBoundaryMesh&>(pointBm).push_back
        (
            new meshPointPatch
            (
                allPointPatchName,
                multiZoneMeshPoints,
                vectorField
                (
                    meshPointNormals,   // is pointNormal expanded to all mesh
                    multiZoneMeshPoints
                ),
                pointBm.size(),
                pointBm,
                meshPointPatch::typeName
            )
        );

        if (verbose)
        {
            const auto& ppp = pointBm.last();
            Info<< "Added multi-patch pointPatch " << ppp.name()
                << " index " << ppp.index()
                << " with " << nMultiPoints << " points"
                << endl;
        }
    }

    const word allEdgePatchName("boundaryEdges");

    const label nPatchEdgePoints =
        returnReduce(twoZoneMeshPoints.size(), sumOp<label>());
    if (nPatchEdgePoints)
    {
        const_cast<pointBoundaryMesh&>(pointBm).push_back
        (
            new meshPointPatch
            (
                allEdgePatchName,
                twoZoneMeshPoints,
                vectorField
                (
                    patchEdgeNormal,    // is pointNormal expanded to all mesh
                    twoZoneMeshPoints
                ),
                pointBm.size(),
                pointBm,
                meshPointPatch::typeName
            )
        );

        if (verbose)
        {
            const auto& ppp = pointBm.last();
            Info<< "Added inter-patch pointPatch " << ppp.name()
                << " index " << ppp.index()
                << " with " << nPatchEdgePoints << " points"
                << endl;
        }
    }

    return pointBm.size() - nPointPatches;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract patch or faceZone surfaces from a polyMesh."
        " The name is historical, it only triangulates faces"
        " when the output format requires it."
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");

    argList::addArgument("output", "The output surface file");

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "excludeProcPatches",
        "Exclude processor patches"
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to extract\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'"
    );
    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to extract.\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "exclude-patches",
        "wordRes",
        "Specify single patch or multiple patches to exclude from -patches."
        " Eg, 'outlet' or '( inlet \".*Wall\" )'",
        true  // mark as an advanced option
    );
    argList::addOptionCompat("exclude-patches", {"excludePatches", 2306});
    argList::addOption
    (
        "featureAngle",
        "angle",
        "Auto-extract feature edges/points and put into separate point-patches"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const auto userOutFileName = args.get<fileName>(1);

    if (!userOutFileName.has_ext())
    {
        FatalErrorInFunction
            << "Missing extension on output name " << userOutFileName
            << exit(FatalError);
    }

    Info<< "Extracting surface from boundaryMesh ..." << nl << nl;

    const bool includeProcPatches =
        (!UPstream::parRun() && !args.found("excludeProcPatches"));

    if (includeProcPatches)
    {
        Info<< "Including all processor patches." << nl << endl;
    }
    else if (UPstream::parRun())
    {
        Info<< "Excluding all processor patches." << nl << endl;
    }

    wordRes includePatches, excludePatches;
    if (args.readListIfPresent<wordRe>("patches", includePatches))
    {
        Info<< "Including patches " << flatOutput(includePatches)
            << nl << endl;
    }
    if (args.readListIfPresent<wordRe>("exclude-patches", excludePatches))
    {
        Info<< "Excluding patches " << flatOutput(excludePatches)
            << nl << endl;
    }

    // Non-mandatory
    const wordRes selectedFaceZones(args.getList<wordRe>("faceZones", false));
    if (selectedFaceZones.size())
    {
        Info<< "Including faceZones " << flatOutput(selectedFaceZones)
            << nl << endl;
    }

    scalar featureAngle = 180.0;
    const bool specifiedFeature = args.readIfPresent
    (
        "featureAngle",
        featureAngle
    );

    Info<< "Reading mesh from time " << runTime.value() << endl;

    #include "createNamedPolyMesh.H"
    if (specifiedFeature)
    {
        Info<< "Detecting all sharp (>" << featureAngle
            << " degrees) patch edges." << nl << endl;

        //#include "createNamedPointMesh.H"
        // Do not read constant/pointMesh - construct from polyMesh only
        Info<< "Create pointMesh for time = "
             << runTime.timeName() << Foam::nl << Foam::endl;
        (void)pointMesh::New(mesh);
    }


    // User specified times
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        Info<< "Time [" << timeIndex << "] = " << runTime.timeName();

        fileName outFileName;
        if (timeDirs.size() == 1)
        {
            outFileName = userOutFileName;
        }
        else
        {
            polyMesh::readUpdateState meshState = mesh.readUpdate();
            if (timeIndex && meshState == polyMesh::UNCHANGED)
            {
                Info<<"  ... no mesh change." << nl;
                continue;
            }

            // The filename based on the original, but with additional
            // time information. The extension was previously checked that
            // it exists
            const auto dot = userOutFileName.rfind('.');

            outFileName =
                userOutFileName.substr(0, dot) + "_"
              + Foam::name(runTime.value()) + "."
              + userOutFileName.ext();
        }

        Info<< nl;

        // Create local surface from:
        // - explicitly named patches only (-patches (at your option)
        // - all patches (default in sequential mode)
        // - all non-processor patches (default in parallel mode)
        // - all non-processor patches (sequential mode, -excludeProcPatches
        //   (at your option)

        // Construct table of patches to include.
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelList patchIds =
        (
            (includePatches.size() || excludePatches.size())
          ? getSelectedPatches(bMesh, includePatches, excludePatches)
          : includeProcPatches
          ? identity(bMesh.size())
          : identity(bMesh.nNonProcessor())
        );

        labelList faceZoneIds;

        const faceZoneMesh& fzm = mesh.faceZones();

        if (selectedFaceZones.size())
        {
            faceZoneIds = fzm.indices(selectedFaceZones);

            Info<< "Additionally extracting faceZones "
                << fzm.names(selectedFaceZones) << nl;
        }


        // From (name of) patch to compact 'zone' index
        HashTable<label> compactZoneID(1024);
        // Mesh face and compact zone indx
        DynamicList<label> faceLabels;
        DynamicList<label> compactZones;
        // Per compact 'zone' index the name and location
        surfZoneIdentifierList surfZones;

        {
            // Collect sizes. Hash on names to handle local-only patches (e.g.
            //  processor patches)
            HashTable<label> patchSize(1024);
            label nFaces = 0;
            for (const label patchi : patchIds)
            {
                const polyPatch& pp = bMesh[patchi];
                patchSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }

            HashTable<label> zoneSize(1024);
            for (const label zonei : faceZoneIds)
            {
                const faceZone& pp = fzm[zonei];
                zoneSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }


            Pstream::mapCombineGather(patchSize, plusEqOp<label>());
            Pstream::mapCombineGather(zoneSize, plusEqOp<label>());

            if (Pstream::master())
            {
                // Allocate compact numbering for all patches/faceZones
                forAllConstIters(patchSize, iter)
                {
                    compactZoneID.insert(iter.key(), compactZoneID.size());
                }

                forAllConstIters(zoneSize, iter)
                {
                    compactZoneID.insert(iter.key(), compactZoneID.size());
                }
            }
            Pstream::broadcast(compactZoneID);


            // Zones
            surfZones.resize_nocopy(compactZoneID.size());
            forAllConstIters(compactZoneID, iter)
            {
                surfZones[*iter] = surfZoneIdentifier(iter.key(), *iter);
                Info<< "surfZone " << *iter
                    <<  " : "      << surfZones[*iter].name()
                    << endl;
            }


            // Rework HashTable into labelList just for speed of conversion
            labelList patchToCompactZone(bMesh.size(), -1);
            labelList faceZoneToCompactZone(bMesh.size(), -1);
            forAllConstIters(compactZoneID, iter)
            {
                label patchi = bMesh.findPatchID(iter.key());
                if (patchi != -1)
                {
                    patchToCompactZone[patchi] = iter.val();
                }
                else
                {
                    label zoneI = fzm.findZoneID(iter.key());
                    faceZoneToCompactZone[zoneI] = iter.val();
                }
            }


            faceLabels.setCapacity(nFaces);
            compactZones.setCapacity(nFaces);

            // Collect faces on patches
            for (const label patchi : patchIds)
            {
                const polyPatch& pp = bMesh[patchi];
                forAll(pp, i)
                {
                    faceLabels.append(pp.start()+i);
                    compactZones.append(patchToCompactZone[pp.index()]);
                }
            }
            // Collect faces on faceZones
            for (const label zonei : faceZoneIds)
            {
                const faceZone& pp = fzm[zonei];
                forAll(pp, i)
                {
                    faceLabels.append(pp[i]);
                    compactZones.append(faceZoneToCompactZone[pp.index()]);
                }
            }
        }


        // Addressing engine for all faces
        const uindirectPrimitivePatch allBoundary
        (
            UIndirectList<face>(mesh.faces(), faceLabels),
            mesh.points()
        );


        // Find correspondence to master points
        labelList pointToGlobal;
        labelList uniqueMeshPoints;
        autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
        (
            allBoundary.meshPoints(),
            allBoundary.meshPointMap(),
            pointToGlobal,
            uniqueMeshPoints
        );

        // Gather all unique points on master
        List<pointField> gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = pointField
        (
            mesh.points(),
            uniqueMeshPoints
        );
        Pstream::gatherList(gatheredPoints);

        // Gather all faces
        List<faceList> gatheredFaces(Pstream::nProcs());
        gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
        for (face& f : gatheredFaces[Pstream::myProcNo()])
        {
            inplaceRenumber(pointToGlobal, f);
        }
        Pstream::gatherList(gatheredFaces);

        // Gather all ZoneIDs
        List<labelList> gatheredZones(Pstream::nProcs());
        gatheredZones[Pstream::myProcNo()] = compactZones;
        Pstream::gatherList(gatheredZones);

        // On master combine all points, faces, zones
        if (Pstream::master())
        {
            pointField allPoints = ListListOps::combine<pointField>
            (
                gatheredPoints,
                accessOp<pointField>()
            );
            gatheredPoints.clear();

            faceList allFaces = ListListOps::combine<faceList>
            (
                gatheredFaces,
                accessOp<faceList>()
            );
            gatheredFaces.clear();

            labelList allZones = ListListOps::combine<labelList>
            (
                gatheredZones,
                accessOp<labelList>()
            );
            gatheredZones.clear();


            UnsortedMeshedSurface<face> unsortedFace
            (
                std::move(allPoints),
                std::move(allFaces),
                std::move(allZones),
                surfZones
            );


            MeshedSurface<face> sortedFace(unsortedFace);

            fileName globalCasePath
            (
                outFileName.isAbsolute()
              ? outFileName
              : (
                    runTime.processorCase()
                  ? runTime.globalPath()/outFileName
                  : runTime.path()/outFileName
                )
            );

            Info<< "Writing merged surface to " << globalCasePath << endl;

            sortedFace.write(globalCasePath);
        }


        if (specifiedFeature)
        {
            // Add edge patches
            const auto& pMesh = pointMesh::New(mesh);

            const label nAdded = addMeshPointPatches
            (
                mesh,
                pMesh,

                allBoundary,    // all patches together
                compactZones,   // originating compactZone
                surfZones,      // per compactZone the index

                featureAngle,
                featureAngle
            );

            if (nAdded)
            {
                pMesh.boundary().write();
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
