/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    surfaceLambdaMuSmooth

Group
    grpSurfaceUtilities

Description
    Smooth a surface using lambda/mu smoothing.

    To get laplacian smoothing, set lambda to the relaxation factor and mu to
    zero.

    Provide an edgeMesh file containing points that are not to be moved during
    smoothing in order to preserve features.

    lambda/mu smoothing: G. Taubin, IBM Research report Rc-19923 (02/01/95)
    "A signal processing approach to fair surface design"

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "boundBox.H"
#include "edgeMesh.H"
#include "matchPoints.H"
#include "MeshedSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<pointField> avg
(
    const meshedSurface& s,
    const bitSet& fixedPoints
)
{
    const labelListList& pointEdges = s.pointEdges();

    auto tavg = tmp<pointField>::New(s.nPoints(), Zero);
    auto& avg = tavg.ref();

    forAll(pointEdges, vertI)
    {
        vector& avgPos = avg[vertI];

        if (fixedPoints.test(vertI))
        {
            avgPos = s.localPoints()[vertI];
        }
        else
        {
            const labelList& pEdges = pointEdges[vertI];

            forAll(pEdges, myEdgeI)
            {
                const edge& e = s.edges()[pEdges[myEdgeI]];

                label otherVertI = e.otherVertex(vertI);

                avgPos += s.localPoints()[otherVertI];
            }

            avgPos /= pEdges.size();
        }
    }

    return tavg;
}


void getFixedPoints
(
    const edgeMesh& feMesh,
    const pointField& points,
    bitSet& fixedPoints
)
{
    scalarList matchDistance(feMesh.points().size(), 1e-1);
    labelList from0To1;

    bool matchedAll = matchPoints
    (
        feMesh.points(),
        points,
        matchDistance,
        false,
        from0To1
    );

    if (!matchedAll)
    {
        WarningInFunction
            << "Did not match all feature points to points on the surface"
            << endl;
    }

    forAll(from0To1, fpI)
    {
        if (from0To1[fpI] != -1)
        {
            fixedPoints.set(from0To1[fpI]);
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Smooth a surface using lambda/mu smoothing.\n"
        "For laplacian smoothing, set lambda to the relaxation factor"
        " and mu to zero."
    );

    argList::noParallel();
    argList::validOptions.clear();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("lambda", "On the interval [0,1]");
    argList::addArgument("mu", "On the interval [0,1]");
    argList::addArgument("iterations", "The number of iterations to perform");
    argList::addArgument("output", "The output surface file");

    argList::addOption
    (
        "featureFile",
        "Fix points from a file containing feature points and edges"
    );
    argList args(argc, argv);

    const auto surfFileName = args.get<fileName>(1);
    const auto lambda = args.get<scalar>(2);
    const auto mu = args.get<scalar>(3);
    const auto iters = args.get<label>(4);
    const auto outFileName = args.get<fileName>(5);

    if (lambda < 0 || lambda > 1)
    {
        FatalErrorInFunction
            << lambda << endl
            << "0: no change   1: move vertices to average of neighbours"
            << exit(FatalError);
    }
    if (mu < 0 || mu > 1)
    {
        FatalErrorInFunction
            << mu << endl
            << "0: no change   1: move vertices to average of neighbours"
            << exit(FatalError);
    }

    Info<< "lambda      : " << lambda << nl
        << "mu          : " << mu << nl
        << "Iters       : " << iters << nl
        << "Reading surface from " << surfFileName << " ..." << endl;

    meshedSurface surf1(surfFileName);

    Info<< "Faces       : " << surf1.size() << nl
        << "Vertices    : " << surf1.nPoints() << nl
        << "Bounding Box: " << boundBox(surf1.localPoints()) << endl;

    bitSet fixedPoints(surf1.localPoints().size(), false);

    if (args.found("featureFile"))
    {
        const auto featureFileName = args.get<fileName>("featureFile");
        Info<< "Reading features from " << featureFileName << " ..." << endl;

        edgeMesh feMesh(featureFileName);

        getFixedPoints(feMesh, surf1.localPoints(), fixedPoints);

        Info<< "Number of fixed points on surface = " << fixedPoints.count()
            << endl;
    }

    for (label iter = 0; iter < iters; iter++)
    {
        // Lambda
        {
            pointField newLocalPoints
            (
                (1 - lambda)*surf1.localPoints()
              + lambda*avg(surf1, fixedPoints)
            );

            pointField newPoints(surf1.points());
            UIndirectList<point>(newPoints, surf1.meshPoints()) =
                newLocalPoints;

            surf1.movePoints(newPoints);
        }

        // Mu
        if (mu != 0)
        {
            pointField newLocalPoints
            (
                (1 + mu)*surf1.localPoints()
              - mu*avg(surf1, fixedPoints)
            );

            pointField newPoints(surf1.points());
            UIndirectList<point>(newPoints, surf1.meshPoints()) =
                newLocalPoints;

            surf1.movePoints(newPoints);
        }
    }

    Info<< "Writing surface to " << outFileName << " ..." << endl;
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
