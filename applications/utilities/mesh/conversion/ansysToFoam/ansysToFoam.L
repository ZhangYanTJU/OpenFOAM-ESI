/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    ansysToFoam

Group
    grpMeshConversionUtilities

Description
    Convert an ANSYS input mesh file (exported from I-DEAS)
    to OpenFOAM format.

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------- *\
   ------ local definitions
 \* ------------------------------------------------------------------------- */

#include <sstream>
// For EOF only
#include <cstdio>

#include "scalar.H"
#include "SpanStream.H"

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "preservePatchTypes.H"
#include "cellShape.H"

// Flex may use register, which is deprecated and incompatible with C++17
#pragma clang diagnostic ignored "-Wdeprecated-register"
// A 'nothing' define to effectively remove from code as well (issue #3337)
#undef  register
#define register

using namespace Foam;

DynamicList<point> dynPoints;
DynamicList<label> dynPointMap;
label maxNodei = 0;

DynamicList<labelList> dynCellLabels;
DynamicList<label> dynCellMap;
DynamicList<label> dynCellType;
label maxCelli = 0;

// Per patch index, the cell/face pairs
PtrList<DynamicList<labelPair>> patchCellFaces;

// Cell types
Map<word> cellTypes;
label currentTypei = -1;

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]

identifier            {alpha}({alpha}|{digit})*
integer               {digit}+
label                 [1-9]{digit}*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}
value                 {floatNum}

node                  ^{space}"N"{cspace}
element               ^{space}"EN"{cspace}
bface                 ^{space}"SFE"{cspace}
elementTypeName       ^{space}"ET"{cspace}
elementType           ^{space}"TYPE"{cspace}


%%

 /* ------------------------------------------------------------------------- *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------- */

{node}{label}{cspace}{x}{cspace}{y}{cspace}{z}{space}\n {
        ISpanStream is(YYText());
        char tag, c;
        label nodei;
        point& node = dynPoints.emplace_back();

        is
            >> tag                  // skip 'N'
            >> c >> nodei
            >> c >> node.x()
            >> c >> node.y()
            >> c >> node.z();

        if (nodei > maxNodei) maxNodei = nodei;

        dynPointMap.push_back(nodei);
    }


{element}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{cspace}{label}{space}\n {
        ISpanStream is(YYText());
        char tag, c;
        label celli;
        labelList& labels = dynCellLabels.emplace_back(8);

        is
            >> tag >> tag           // skip 'EN'
            >> c >> celli
            >> c >> labels[0]
            >> c >> labels[1]
            >> c >> labels[2]
            >> c >> labels[3]
            >> c >> labels[4]
            >> c >> labels[5]
            >> c >> labels[6]
            >> c >> labels[7];

        if (celli > maxCelli) maxCelli = celli;

        dynCellMap.push_back(celli);
        dynCellType.push_back(currentTypei);
    }


{bface}{label}{cspace}{label}{cspace}{identifier}{cspace}{integer}{cspace}{value}{space}\n {
        ISpanStream is(YYText());
        char tag, c;
        label elementi, facei;
        scalar indexValue, unknown;

        is
            >> tag >> tag >> tag    // skip 'SFE'
            >> c >> elementi
            >> c >> facei
            >> c >> tag >> tag >> tag >> tag
            >> c >> unknown
            >> c >> indexValue;

        const label patchi = label(indexValue);

        if (patchCellFaces.size() < patchi)
        {
            label i = patchCellFaces.size();

            patchCellFaces.resize(patchi);

            for (; i < patchi; ++i)
            {
                patchCellFaces.try_emplace(i);
            }
        }

        patchCellFaces[patchi-1].emplace_back(elementi, facei);
    }


{elementTypeName}{label}{cspace}{identifier}{space}\n {
        ISpanStream is(YYText());
        char tag,c;
        label cellTypei;
        word cellTypeName;

        is
            >> tag >> tag           // skip 'ET'
            >> c >> cellTypei
            >> c >> cellTypeName;

        Info<< "Read typeName " << cellTypeName
            << " for type " << cellTypei << endl;

        cellTypes.insert(cellTypei, cellTypeName);
    }


{elementType}{label}{space}\n {
        ISpanStream is(YYText());
        char tag,c;
        label cellTypei;

        is
            >> tag >> tag >> tag >> tag     // skip 'TYPE'
            >> c >> cellTypei;

        currentTypei = cellTypei;
    }



 /* ------------------------------------------------------------------------- *\
    ------ Ignore remaining space and \n s.  Any other characters are errors.
 \* ------------------------------------------------------------------------- */

.|\n {}


 /* ------------------------------------------------------------------------- *\
    ------ On EOF return to previous file, if none exists terminate.
 \* ------------------------------------------------------------------------- */

<<EOF>> {
            yyterminate();
    }
%%


#include "fileName.H"
#include <fstream>


label findFace(const polyMesh& mesh, const face& f)
{
    const labelList& pFaces = mesh.pointFaces()[f[0]];

    for (const label facei : pFaces)
    {
        if (mesh.faces()[facei] == f)
        {
            return facei;
        }
    }

    FatalErrorInFunction
        << "Cannot find a face matching " << f
        << exit(FatalError);

    return -1;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert an ANSYS input mesh file (exported from I-DEAS)"
        " to OpenFOAM format"
    );

    argList::noParallel();
    argList::addArgument("ANSYS input file");
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1"
    );

    argList args(argc, argv);

    if (!args.check())
    {
        FatalError.exit();
    }

    // Actually uses default=0 to skip unnecessary scaling by factor 1.
    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 0);

    #include "createTime.H"

    const auto ansysFile = args.get<fileName>(1);
    std::ifstream ansysStream(ansysFile);

    if (!ansysStream)
    {
        FatalErrorInFunction
            << args.executable()
            << ": file " << ansysFile << " not found"
            << exit(FatalError);
    }

    yyFlexLexer lexer(&ansysStream);
    while (lexer.yylex() != 0)
    {}

    Info<< "Creating points" << endl;

    pointField points(std::move(dynPoints));

    // Scale points by the given scale factor
    if (scaleFactor > 0)
    {
        points *= scaleFactor;
    }

    labelList pointMap(maxNodei+1);
    {
        label i = 0;
        for (const label pointi : dynPointMap)
        {
            pointMap[pointi] = i++;
        }
    }

    Info<< "Creating cells" << endl;

    labelList cellMap(maxCelli+1);
    {
        label i = 0;
        for (const label celli : dynCellMap)
        {
            cellMap[celli] = i++;
        }
    }


    const cellModel& hex = cellModel::ref(cellModel::HEX);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& pyr = cellModel::ref(cellModel::PYR);
    const cellModel& tet = cellModel::ref(cellModel::TET);

    FixedList<label, 8> labelsHex;
    FixedList<label, 6> labelsPrism;
    FixedList<label, 5> labelsPyramid;
    FixedList<label, 4> labelsTet;

    cellShapeList cellShapes(dynCellLabels.size());
    label nCells = 0;

    for (const labelList& labels : dynCellLabels)
    {
        if      // Tetrahedron
        (
            labels[2] == labels[3]
         && labels[4] == labels[5]
         && labels[5] == labels[6]
         && labels[6] == labels[7]
        )
        {
            labelsTet[0] = pointMap[labels[0]];
            labelsTet[1] = pointMap[labels[1]];
            labelsTet[2] = pointMap[labels[2]];
            labelsTet[3] = pointMap[labels[4]];

            cellShapes[nCells++].reset(tet, labelsTet);
        }

        else if // Square-based pyramid
        (
            labels[4] == labels[5]
         && labels[5] == labels[6]
         && labels[6] == labels[7]
        )
        {
            labelsPyramid[0] = pointMap[labels[0]];
            labelsPyramid[1] = pointMap[labels[1]];
            labelsPyramid[2] = pointMap[labels[2]];
            labelsPyramid[3] = pointMap[labels[3]];
            labelsPyramid[4] = pointMap[labels[4]];

            cellShapes[nCells++].reset(pyr, labelsPyramid);
        }

        else if // Triangular prism
        (
            labels[2] == labels[3]
         && labels[6] == labels[7]
        )
        {
            labelsPrism[0] = pointMap[labels[0]];
            labelsPrism[1] = pointMap[labels[1]];
            labelsPrism[2] = pointMap[labels[2]];
            labelsPrism[3] = pointMap[labels[4]];
            labelsPrism[4] = pointMap[labels[5]];
            labelsPrism[5] = pointMap[labels[6]];

            cellShapes[nCells++].reset(prism, labelsPrism);
        }

        else // Hex
        {
            labelsHex[0] = pointMap[labels[0]];
            labelsHex[1] = pointMap[labels[1]];
            labelsHex[2] = pointMap[labels[2]];
            labelsHex[3] = pointMap[labels[3]];
            labelsHex[4] = pointMap[labels[4]];
            labelsHex[5] = pointMap[labels[5]];
            labelsHex[6] = pointMap[labels[6]];
            labelsHex[7] = pointMap[labels[7]];

            cellShapes[nCells++].reset(hex, labelsHex);
        }
    }


    const word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;

    // Create dummy mesh just to find out what are internal/external faces
    auto dummyMesh = autoPtr<polyMesh>::New
    (
        IOobject
        (
            "dummyMesh",
            runTime.constant(),
            runTime
        ),
        pointField(points),  // copy
        cellShapes,
        faceListList(),
        wordList(),
        wordList(),
        defaultFacesName,
        defaultFacesType,
        wordList()
    );


    // Warning: tet face order has changed between version 1.9.6 and 2.0
    //
    const label faceIndex[7][6] =
    {
        {-1, -1, -1, -1, -1, -1}, // 0
        {-1, -1, -1, -1, -1, -1}, // 1
        {-1, -1, -1, -1, -1, -1}, // 2
        {-1, -1, -1, -1, -1, -1}, // 3
        { 3,  2,  0, -1,  1, -1}, // tet (version 2.0)
        { 0,  4,  3, -1,  2,  1}, // prism
        { 4,  2,  1,  3,  0,  5}, // hex
    };

    Info<< "Creating boundary patches" << endl;

    faceListList boundary(patchCellFaces.size());
    wordList patchNames(patchCellFaces.size());

    forAll(patchCellFaces, patchi)
    {
        DynamicList<face> patchFaces;

        patchFaces.reserve(patchCellFaces[patchi].size());

        for (const auto& tup : patchCellFaces[patchi])
        {
            const label celli = tup.first();
            const label facei = tup.second();

            const cellShape& shape = cellShapes[cellMap[celli]];

            patchFaces.push_back
            (
                shape.faces()[faceIndex[shape.nFaces()][facei-1]]
            );
        }

        boundary[patchi] = std::move(patchFaces);
        patchNames[patchi] = polyPatch::defaultName(patchi + 1);
    }


    //
    // Lookup the face labels for all the boundary faces
    //
    labelListList boundaryFaceLabels(boundary.size());
    forAll(boundary, patchi)
    {
        const faceList& bFaces = boundary[patchi];
        labelList& bFaceLabels = boundaryFaceLabels[patchi];
        bFaceLabels.setSize(bFaces.size());
        forAll(bFaces, i)
        {
            bFaceLabels[i] = findFace(dummyMesh(), bFaces[i]);
        }
    }


    // Now split the boundary faces into external and internal faces. All
    // faces go into faceZones and external faces go into patches.
    List<faceList> patchFaces(patchCellFaces.size());
    labelList patchNFaces(patchCellFaces.size(), Zero);
    forAll(boundary, patchi)
    {
        const faceList& bFaces = boundary[patchi];
        const labelList& bFaceLabels = boundaryFaceLabels[patchi];

        patchFaces[patchi].setSize(bFaces.size());

        forAll(bFaces, i)
        {
            if (!dummyMesh().isInternalFace(bFaceLabels[i]))
            {
                patchFaces[patchi][patchNFaces[patchi]++] = bFaces[i];
            }
        }
        patchFaces[patchi].setSize(patchNFaces[patchi]);

        Info<< "Patch " << patchi << " named " << patchNames[patchi]
            << ": " << boundary[patchi].size() << " faces" << endl;
    }

    // We no longer need the dummyMesh
    dummyMesh.clear();


    Info<< "ansysToFoam: " << endl
        << "Ansys file format does not provide information about the type of "
        << "the patch (eg. wall, symmetry plane, cyclic etc)." << endl
        << "All the patches have been created "
        << "as type patch. Please reset after mesh conversion as necessary."
        << endl;

    PtrList<dictionary> patchDicts;

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );

    // Add information to dictionary
    forAll(patchNames, patchi)
    {
        if (!patchDicts.set(patchi))
        {
            patchDicts.set(patchi, new dictionary());
        }
        // Add but not overwrite
        patchDicts[patchi].add("type", polyPatch::typeName, false);
    }


    polyMesh pShapeMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points),
        cellShapes,
        patchFaces,
        patchNames,
        patchDicts,
        defaultFacesName,
        defaultFacesType
    );


    if (cellTypes.size() || patchNames.size())
    {
        DynamicList<pointZone*> pz;
        DynamicList<faceZone*> fz;
        DynamicList<cellZone*> cz;

        // FaceZones
        forAll(boundaryFaceLabels, patchi)
        {
            if (boundaryFaceLabels[patchi].size())
            {
                // Re-do the boundaryFaceLabels since the boundary face
                // labels will be different on the pShapeMesh.
                const faceList& bFaces = boundary[patchi];
                labelList& bFaceLabels = boundaryFaceLabels[patchi];
                forAll(bFaceLabels, i)
                {
                    bFaceLabels[i] = findFace(pShapeMesh, bFaces[i]);
                }

                Info<< "Creating faceZone " <<  patchNames[patchi]
                    << " with " << bFaceLabels.size() << " faces" << endl;

                fz.push_back
                (
                    new faceZone
                    (
                        patchNames[patchi],
                        bFaceLabels,
                        false, // none are flipped
                        fz.size(),
                        pShapeMesh.faceZones()
                    )
                );
            }
        }


        // CellZones
        for (const label cellType : cellTypes.sortedToc())
        {
            // Pick up cells in zone
            DynamicList<label> addr;

            auto cellMapIter = dynCellMap.cbegin();

            for (const auto& ctype : dynCellType)
            {
                if (ctype == cellType)
                {
                    addr.push_back(cellMap[*cellMapIter]);
                }
                ++cellMapIter;
            }

            Info<< "Creating cellZone " << cellTypes[cellType]
                << " with " << addr.size() << " cells" << endl;

            cz.append
            (
                new cellZone
                (
                    cellTypes[cellType],
                    addr,
                    cz.size(),
                    pShapeMesh.cellZones()
                )
            );
        }

        pShapeMesh.addZones(pz, fz, cz);
    }


    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing polyMesh" << endl;
    pShapeMesh.removeFiles();
    pShapeMesh.write();

    Info<< nl << "End" << endl;
    return 0;
}


 /* ------------------------------------------------------------------------- *\
    ------ End of ansysToFoam.L
 \* ------------------------------------------------------------------------- */
