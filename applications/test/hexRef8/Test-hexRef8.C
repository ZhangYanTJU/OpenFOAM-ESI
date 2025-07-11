/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Test-hexRef8

Description
    Test app for refinement and unrefinement. Runs a few iterations refining
    and unrefining.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "Random.H"
#include "calculatedPointPatchFields.H"
#include "pointConstraints.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addBoolOption
    (
        "inflate",
        "Use inflation/deflation for splitting/deleting cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // Allow override of time from specified time options, or no-op
    timeSelector::setTimeIfPresent(runTime, args);

    #include "createMesh.H"

    const bool inflate = args.found("inflate");

    if (inflate)
    {
        Info<< "Splitting/deleting cells using inflation/deflation"
            << nl << endl;
    }
    else
    {
        Info<< "Splitting/deleting cells, introducing points at new position"
            << nl << endl;
    }


    const pointConstraints& pc = pointConstraints::New(pointMesh::New(mesh));

    Random rndGen(0);


    // Force generation of V()
    (void)mesh.V();


    // Test mapping
    // ------------

    // 1. uniform field stays uniform
    volScalarField one
    (
        IOobject
        (
            "one",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0),
        fvPatchFieldBase::zeroGradientType()
    );
    Info<< "Writing one field "
        << one.name() << " in " << runTime.timeName() << endl;
    one.write();


    // 2. linear profile gets preserved
    volScalarField ccX
    (
        IOobject
        (
            "ccX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C().component(0)
    );
    Info<< "Writing x component of cell centres to "
        << ccX.name()
        << " in " << runTime.timeName() << endl;
    ccX.write();


    // Uniform surface field
    surfaceScalarField surfaceOne
    (
        IOobject
        (
            "surfaceOne",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );
    Info<< "Writing surface one field "
        << surfaceOne.name() << " in " << runTime.timeName() << endl;
    surfaceOne.write();


    // Uniform point field
    pointScalarField pointX
    (
        IOobject
        (
            "pointX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar("one", dimless, 1.0)
    );
    pointX.primitiveFieldRef() = mesh.points().component(0);
    pointX.correctBoundaryConditions();
    Info<< "Writing x-component field "
        << pointX.name() << " in " << runTime.timeName() << endl;
    pointX.write();



    // Force allocation of V. Important for any mesh changes since otherwise
    // old time volumes are not stored
    const scalar totalVol = gSum(mesh.V());


    // Construct refiner. Read initial cell and point levels.
    hexRef8 meshCutter(mesh);

    // Comparison for inequality
    const auto isNotEqual = notEqualOp<scalar>(1e-10);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (mesh.globalData().nTotalCells() == 0)
        {
            break;
        }


        mesh.moving(false);
        mesh.topoChanging(false);

        label action = rndGen.position<label>(0, 5);

        if (action == 0)
        {
            Info<< nl << "-- moving only" << endl;
            mesh.movePoints(pointField(mesh.points()));
        }
        else if (action == 1 || action == 2)
        {
            // Mesh changing engine.
            polyTopoChange meshMod(mesh);

            if (action == 1)
            {
                // Refine
                label nRefine = mesh.nCells()/20;
                DynamicList<label> refineCandidates(nRefine);

                for (label i=0; i<nRefine; i++)
                {
                    refineCandidates.append
                    (
                        rndGen.position<label>(0, mesh.nCells()-1)
                    );
                }

                labelList cellsToRefine
                (
                    meshCutter.consistentRefinement
                    (
                        refineCandidates,
                        true                  // buffer layer
                    )
                );
                Info<< nl << "-- selected "
                    << returnReduce(cellsToRefine.size(), sumOp<label>())
                    << " cells out of " << mesh.globalData().nTotalCells()
                    << " for refinement" << endl;

                // Play refinement commands into mesh changer.
                meshCutter.setRefinement(cellsToRefine, meshMod);
            }
            else
            {
                // Unrefine
                labelList allSplitPoints(meshCutter.getSplitPoints());

                label nUnrefine = allSplitPoints.size()/20;
                labelHashSet candidates(2*nUnrefine);

                for (label i=0; i<nUnrefine; i++)
                {
                    const label index =
                        rndGen.position<label>(0, allSplitPoints.size()-1);

                    candidates.insert(allSplitPoints[index]);
                }

                labelList splitPoints = meshCutter.consistentUnrefinement
                (
                    candidates.toc(),
                    false
                );
                Info<< nl << "-- selected "
                    << returnReduce(splitPoints.size(), sumOp<label>())
                    << " points out of "
                    << returnReduce(allSplitPoints.size(), sumOp<label>())
                    << " for unrefinement" << endl;

                // Play refinement commands into mesh changer.
                meshCutter.setUnrefinement(splitPoints, meshMod);
            }




            // Create mesh, return map from old to new mesh.
            Info<< nl << "-- actually changing mesh" << endl;
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, inflate);

            // Update fields
            Info<< nl << "-- mapping mesh data" << endl;
            mesh.updateMesh(map());

            // Inflate mesh
            if (map().hasMotionPoints())
            {
                Info<< nl << "-- moving mesh" << endl;
                mesh.movePoints(map().preMotionPoints());
            }

            // Update numbering of cells/vertices.
            Info<< nl << "-- mapping hexRef8 data" << endl;
            meshCutter.updateMesh(map());
        }


        Info<< nl<< "-- Mesh : moving:" << mesh.moving()
            << " topoChanging:" << mesh.topoChanging()
            << " changing:" << mesh.changing()
            << endl;



        Info<< "Writing fields" << nl << endl;
        runTime.write();



        // Check mesh volume conservation
        if (mesh.moving())
        {
            #include "volContinuity.H"
        }
        else
        {
            if (mesh.V().size() != mesh.nCells())
            {
                FatalErrorInFunction
                    << "Volume not mapped. V:" << mesh.V().size()
                    << " nCells:" << mesh.nCells()
                    << exit(FatalError);
            }

            const scalar newVol = gSum(mesh.V());
            Info<< "Initial volume = " << totalVol
                << "  New volume = " << newVol
                << endl;

            if (mag(newVol-totalVol)/totalVol > 1e-10)
            {
                FatalErrorInFunction
                    << "Volume loss: old volume:" << totalVol
                    << "  new volume:" << newVol
                    << exit(FatalError);
            }
            else
            {
                Info<< "Volume check OK" << nl << endl;
            }
        }


        // Check constant profile
        {
            auto limits =  gMinMax(one);

            Info<< "Uniform one field min = "
                << limits.min() << "  max = " << limits.max() << endl;

            if (isNotEqual(limits.min(), 1) || isNotEqual(limits.max(), 1))
            {
                FatalErrorInFunction
                    << "Uniform volVectorField not preserved."
                    << " Min and max should both be 1.0. min:" << min
                    << " max:" << max
                    << exit(FatalError);
            }
            else
            {
                Info<< "Uniform field mapping check OK" << nl << endl;
            }
        }

        // Check linear profile
        {
            const scalarField diff = ccX-mesh.C().component(0);

            auto limits =  gMinMax(diff);

            Info<< "Linear profile field min = "
                << limits.min() << "  max = " << limits.max() << endl;

            if (isNotEqual(limits.min(), 0) || isNotEqual(limits.max(), 0))
            {
                Info<< "Linear profile not preserved."
                    << " Min and max should both be 0.0. min:" << min
                    << " max:" << max << nl << endl;
            }
            else
            {
                Info<< "Linear profile mapping check OK" << nl << endl;
            }
        }

        // Check face field mapping
        if (surfaceOne.size())
        {
            auto limits =  gMinMax(surfaceOne.primitiveField());

            Info<< "Uniform surface field min = "
                << limits.min() << "  max = " << limits.max() << endl;

            if (isNotEqual(limits.min(), 1) || isNotEqual(limits.max(), 1))
            {
                FatalErrorInFunction
                    << "Uniform surfaceScalarField not preserved."
                    << " Min and max should both be 1.0. min:" << min
                    << " max:" << max
                    << exit(FatalError);
            }
            else
            {
                Info<< "Uniform surfaceScalarField mapping check OK" << nl
                    << endl;
            }
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "pc:" << pc.patchPatchPointConstraintPoints().size() << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
