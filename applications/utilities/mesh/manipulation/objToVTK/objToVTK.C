/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    objToVTK

Group
    grpMeshManipulationUtilities

Description
    Read obj line (not surface) file and convert into legacy VTK file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "stringOps.H"
#include "point.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string getLine(std::ifstream& is)
{
    string line;
    do
    {
        std::getline(is, line);
    }
    while (line.starts_with('#'));

    return line;
}


// Token list with one of the following:
//     f v1 v2 v3 ...
//     f v1/vt1 v2/vt2 v3/vt3 ...
//     l v1 v2 v3 ...
//     l v1/vt1 v2/vt2 v3/vt3 ...
static label readObjVertices
(
    const SubStrings& tokens,
    DynamicList<label>& verts
)
{
    verts.clear();

    bool first = true;
    for (const auto& tok : tokens)
    {
        if (first)
        {
            // skip initial "f" or "l"
            first = false;
            continue;
        }

        std::string vrtSpec(tok.str());

        if
        (
            const auto slash = vrtSpec.find('/');
            slash != std::string::npos
        )
        {
            vrtSpec.erase(slash);
        }

        label vertId = readLabel(vrtSpec);

        verts.push_back(vertId - 1);
    }

    return verts.size();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read obj line (not surface) file and convert into legacy VTK file"
    );

    argList::noParallel();
    argList::addArgument("obj-file", "The input obj line file");
    argList::addArgument("vtk-file", "The output vtk file");
    argList args(argc, argv);

    const auto objName = args.get<fileName>(1);
    const auto outName = args.get<fileName>(2);

    std::ifstream OBJfile(objName);

    Info<< "Processing file " << objName << endl;

    if (!OBJfile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << objName << exit(FatalError);
    }

    // Points and lines
    DynamicList<point> points;
    DynamicList<vector> pointNormals;
    DynamicList<labelList> polyLines;
    DynamicList<labelList> polygons;

    DynamicList<label> dynVerts;

    bool hasWarned = false;

    label lineNo = 0;
    while (OBJfile.good())
    {
        const string line = getLine(OBJfile);
        lineNo++;

        if (line.empty()) continue;

        const auto tokens = stringOps::splitSpace(line);

        // Require command and some arguments
        if (tokens.size() < 2)
        {
            continue;
        }

        const word cmd = word::validate(tokens[0]);

        if (cmd == "v")
        {
            // Vertex
            // v x y z

            points.emplace_back
            (
                readScalar(tokens[1]),
                readScalar(tokens[2]),
                readScalar(tokens[3])
            );
        }
        else if (cmd == "vn")
        {
            // Vertex normals
            // vn x y z

            pointNormals.emplace_back
            (
                readScalar(tokens[1]),
                readScalar(tokens[2]),
                readScalar(tokens[3])
            );
        }
        else if (cmd == "l")
        {
            // Line
            // l v1 v2 v3 ...
            // OR
            // l v1/vt1 v2/vt2 v3/vt3 ...

            readObjVertices(tokens, dynVerts);
            polyLines.emplace_back() = dynVerts;
        }
        else if (cmd == "f")
        {
            // f v1 v2 v3 ...
            // OR
            // f v1/vt1 v2/vt2 v3/vt3 ...

            readObjVertices(tokens, dynVerts);
            polygons.emplace_back() = dynVerts;
        }
        else if (cmd != "")
        {
            if (!hasWarned)
            {
                hasWarned = true;

                WarningInFunction
                    << "Unrecognized OBJ command " << cmd << nl
                    << "In line " << line
                    << " at linenumber " << lineNo << nl
                    << "Only recognized commands are 'v' and 'l'.\n"
                    << "If this is a surface command use surfaceConvert instead"
                    << " to convert to a file format that can be read by VTK"
                    << endl;
            }
        }
    }


    //
    // Write as vtk 'polydata' file
    //


    OFstream outFile(outName);

    outFile
        << "# vtk DataFile Version 2.0\n"
        << objName << nl
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << points.size() << " double\n";

    for (const point& pt : points)
    {
        outFile
            << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }

    outFile
        << "VERTICES " << points.size() << ' ' << (2 * points.size()) << nl;

    forAll(points, i)
    {
        outFile << 1 << ' ' << i << nl;
    }

    label nItems = 0;
    for (const labelList& line : polyLines)
    {
        nItems += line.size() + 1;
    }

    outFile
        << "LINES " << polyLines.size() << ' ' << nItems << nl;

    for (const labelList& line : polyLines)
    {
        outFile << line.size();

        for (const label vrt : line)
        {
            outFile << ' ' << vrt;
        }
        outFile << nl;
    }


    nItems = 0;
    for (const labelList& line : polygons)
    {
        nItems += line.size() + 1;
    }

    outFile
        << "POLYGONS " << polygons.size() << ' ' << nItems << nl;

    for (const labelList& line : polygons)
    {
        outFile << line.size();

        for (const label vrt : line)
        {
            outFile << ' ' << vrt;
        }
        outFile << nl;
    }


    outFile
        << "POINT_DATA " << points.size() << nl
        << "SCALARS pointID double 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(points, i)
    {
        outFile << i;

        if ((i % 10) == 1)
        {
            outFile << nl;
        }
        else
        {
            outFile << ' ';
        }
    }

    if (!pointNormals.empty())
    {
        outFile << nl << "NORMALS pointNormals double\n";

        for(const vector& n : pointNormals)
        {
            outFile
                << float(n.x()) << ' '
                << float(n.y()) << ' '
                << float(n.z()) << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
