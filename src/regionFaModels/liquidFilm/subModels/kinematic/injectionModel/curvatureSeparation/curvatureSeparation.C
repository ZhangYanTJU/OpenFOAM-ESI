/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "curvatureSeparation.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

#include "stringListOps.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(curvatureSeparation, 0);
addToRunTimeSelectionTable
(
    injectionModel,
    curvatureSeparation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<areaScalarField> curvatureSeparation::calcInvR1
(
    const areaVectorField& U
) const
{
    // method 1
/*
    tmp<areaScalarField> tinvR1
    (
        new areaScalarField("invR1", fvc::div(film().nHat()))
    );
*/

    // method 2
    dimensionedScalar smallU("smallU", dimVelocity, ROOTVSMALL);
    areaVectorField UHat(U/(mag(U) + smallU));
    tmp<areaScalarField> tinvR1
    (
        new areaScalarField("invR1", UHat & (UHat & gradNHat_))
    );

    scalarField& invR1 = tinvR1.ref().primitiveFieldRef();

    // apply defined patch radii
    const scalar rMin = 1e-6;
    scalar definedInvR1 = 1.0/max(rMin, definedPatchRadii_);

    //UIndirectList<scalar>(invR1, pbm[patchi].faceCells()) = definedInvR1;

    if (definedPatchRadii_ > 0)
    {
        invR1 = definedInvR1;
    }

    // filter out large radii
    const scalar rMax = 1e6;
    forAll(invR1, i)
    {
        if (mag(invR1[i]) < 1/rMax)
        {
            invR1[i] = -1.0;
        }
    }

    const faMesh& mesh = film().regionMesh();

    if (debug && mesh.time().writeTime())
    {
        tinvR1().write();
    }

    return tinvR1;
}


tmp<scalarField> curvatureSeparation::calcCosAngle
(
    const edgeScalarField& phi
) const
{
    const faMesh& mesh = film().regionMesh();
    const areaVectorField nf(mesh.faceAreaNormals());
    const labelUList& own = mesh.edgeOwner();
    const labelUList& nbr = mesh.edgeNeighbour();

    scalarField phiMax(mesh.nFaces(), -GREAT);
    scalarField cosAngle(mesh.nFaces(), Zero);
    forAll(nbr, edgei)
    {
        label cellO = own[edgei];
        label cellN = nbr[edgei];

        if (phi[edgei] > phiMax[cellO])
        {
            phiMax[cellO] = phi[edgei];
            cosAngle[cellO] = -gHat_ & nf[edgei];
        }
        if (-phi[edgei] > phiMax[cellN])
        {
            phiMax[cellN] = -phi[edgei];
            cosAngle[cellN] = -gHat_ & -nf[edgei];
        }
    }

    forAll(phi.boundaryField(), patchi)
    {
        const faePatchScalarField& phip = phi.boundaryField()[patchi];
        const faPatch& pp = phip.patch();
        const labelList& edgeFaces = pp.edgeFaces();
        const vectorField nf(pp.edgeNormals());
        forAll(phip, i)
        {
            label facei = edgeFaces[i];
            if (phip[i] > phiMax[facei])
            {
                phiMax[facei] = phip[i];
                cosAngle[facei] = -gHat_ & nf[i];
            }
        }
    }
/*
    // correction for cyclics - use cyclic pairs' face normal instead of
    // local face normal
    const fvBoundaryMesh& pbm = mesh.boundary();
    forAll(phi.boundaryField(), patchi)
    {
        if (isA<cyclicPolyPatch>(pbm[patchi]))
        {
            const scalarField& phip = phi.boundaryField()[patchi];
            const vectorField nf(pbm[patchi].nf());
            const labelList& faceCells = pbm[patchi].faceCells();
            const label sizeBy2 = pbm[patchi].size()/2;

            for (label face0=0; face0<sizeBy2; face0++)
            {
                label face1 = face0 + sizeBy2;
                label cell0 = faceCells[face0];
                label cell1 = faceCells[face1];

                // flux leaving half 0, entering half 1
                if (phip[face0] > phiMax[cell0])
                {
                    phiMax[cell0] = phip[face0];
                    cosAngle[cell0] = -gHat_ & -nf[face1];
                }

                // flux leaving half 1, entering half 0
                if (-phip[face1] > phiMax[cell1])
                {
                    phiMax[cell1] = -phip[face1];
                    cosAngle[cell1] = -gHat_ & nf[face0];
                }
            }
        }
    }
*/
    // checks
    if (debug && mesh.time().writeTime())
    {
        areaScalarField volCosAngle
        (
            IOobject
            (
                "cosAngle",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            film().regionMesh(),
            dimensionedScalar(dimless, Zero)
        );
        volCosAngle.primitiveFieldRef() = cosAngle;
        volCosAngle.correctBoundaryConditions();
        volCosAngle.write();
    }

    return max(min(cosAngle, scalar(1)), scalar(-1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureSeparation::curvatureSeparation
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    gradNHat_(fac::grad(film.regionMesh().faceAreaNormals())),
    deltaByR1Min_(coeffDict_.getOrDefault<scalar>("deltaByR1Min", 0)),
    definedPatchRadii_(coeffDict_.getOrDefault<scalar>("definedPatchRadii", 0)),
    magG_(mag(film.g().value())),
    gHat_(Zero)
{
    if (magG_ < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Acceleration due to gravity must be non-zero"
            << exit(FatalError);
    }

    gHat_ = film.g().value()/magG_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

curvatureSeparation::~curvatureSeparation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void curvatureSeparation::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    //const kinematicSingleLayer& film =
    //    refCast<const kinematicSingleLayer>(this->film());

    const faMesh& mesh = film().regionMesh();

    const areaScalarField& delta = film().h();
    const areaVectorField& U = film().Uf();
    const edgeScalarField& phi = film().phi2s();
    const areaScalarField& rho = film().rho();
    const scalarField magSqrU(magSqr(film().Uf()));
    const areaScalarField& sigma = film().sigma();

    const scalarField invR1(calcInvR1(U));
    const scalarField cosAngle(calcCosAngle(phi));

    // calculate force balance
    const scalar Fthreshold = 1e-10;
    scalarField Fnet(mesh.nFaces(), Zero);
    scalarField separated(mesh.nFaces(), Zero);
    forAll(invR1, i)
    {
        if ((invR1[i] > 0) && (delta[i]*invR1[i] > deltaByR1Min_))
        {
            scalar R1 = 1.0/(invR1[i] + ROOTVSMALL);
            scalar R2 = R1 + delta[i];

            // inertial force
            scalar Fi = -delta[i]*rho[i]*magSqrU[i]*72.0/60.0*invR1[i];

            // body force
            scalar Fb =
              - 0.5*rho[i]*magG_*invR1[i]*(sqr(R1) - sqr(R2))*cosAngle[i];

            // surface force
            scalar Fs = sigma[i]/R2;

            Fnet[i] = Fi + Fb + Fs;

            if (Fnet[i] + Fthreshold < 0)
            {
                separated[i] = 1.0;
            }
        }
    }

    // inject all available mass
    massToInject = separated*availableMass;
    diameterToInject = separated*delta;
    availableMass -= separated*availableMass;

    addToInjectedMass(sum(separated*availableMass));

    if (debug && mesh.time().writeTime())
    {
        areaScalarField volFnet
        (
            IOobject
            (
                "Fnet",
                film().primaryMesh().time().timeName(),
                film().primaryMesh(),
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimForce, Zero)
        );
        volFnet.primitiveFieldRef() = Fnet;
        volFnet.write();
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
