/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "componentDisplacementMotionSolver.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(componentDisplacementMotionSolver, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::direction Foam::componentDisplacementMotionSolver::cmpt
(
    const word& cmptName
) const
{
    if (cmptName == "x")
    {
        return vector::X;
    }
    else if (cmptName == "y")
    {
        return vector::Y;
    }
    else if (cmptName == "z")
    {
        return vector::Z;
    }
    else
    {
        FatalErrorInFunction
            << "Given component name " << cmptName << " should be x, y or z"
            << exit(FatalError);

        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::componentDisplacementMotionSolver::componentDisplacementMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const word& type
)
:
    motionSolver(mesh, dict, type),
    cmptName_(coeffDict().get<word>("component")),
    cmpt_(cmpt(cmptName_)),
    points0_
    (
        pointIOField
        (
            IOobject
            (
                "points",
                time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).component(cmpt_)
    ),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement" + cmptName_,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    )
{
    if (points0_.size() != mesh.nPoints())
    {
        const fileName fName
        (
            IOobject
            (
                "points",
                mesh.time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ).typeFilePath<pointIOField>()
        );

        FatalErrorInFunction
            << "Number of points in mesh " << mesh.nPoints()
            << " differs from number of points " << points0_.size()
            << " read from file " << fName << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::componentDisplacementMotionSolver::~componentDisplacementMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::componentDisplacementMotionSolver::movePoints(const pointField& p)
{
    // No local data to update
}


void Foam::componentDisplacementMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    // pointMesh already updates pointFields.

    motionSolver::updateMesh(mpm);

    // Map points0_. Bit special since we somehow have to come up with
    // a sensible points0 position for introduced points.
    // Find out scaling between points0 and current points

    // Get the new points either from the map or the mesh
    const scalarField points
    (
        mpm.hasMotionPoints()
      ? mpm.preMotionPoints().component(cmpt_)
      : mesh().points().component(cmpt_)
    );

    // Get extents of points0 and points and determine scale
    const scalar scale = gMinMax(points0_).span() / gMinMax(points).span();

    scalarField newPoints0(mpm.pointMap().size());

    forAll(newPoints0, pointi)
    {
        label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            label masterPointi = mpm.reversePointMap()[oldPointi];

            if (masterPointi == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
            else
            {
                // New point. Assume motion is scaling.
                newPoints0[pointi] =
                    points0_[oldPointi]
                  + scale*(points[pointi]-points[masterPointi]);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot work out coordinates of introduced vertices."
                << " New vertex " << pointi << " at coordinate "
                << points[pointi] << exit(FatalError);
        }
    }
    points0_.transfer(newPoints0);
}


// ************************************************************************* //
