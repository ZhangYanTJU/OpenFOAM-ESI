/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "PDRkEpsilon.H"
#include "PDRDragModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PDRkEpsilon, 0);
addToRunTimeSelectionTable(RASModel, PDRkEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PDRkEpsilon::PDRkEpsilon
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const fluidThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    Foam::RASModels::kEpsilon<EddyDiffusivity<compressible::turbulenceModel>>
    (
        geometricOneField(),
        rho,
        U,
        phi,
        phi,
        thermophysicalModel,
        turbulenceModelName,
        modelName
    ),

    C5_(coeffDict_.get<scalar>("C5")),
    C6_(coeffDict_.get<scalar>("C6")),
    maxLOverCellW_(coeffDict_.get<scalar>("maxLOverCellW")),
    lCoef_(coeffDict_.get<scalar>("lCoef")),
    noTurbUntil_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "noTurbUntil",
            coeffDict_,
            0.0
        )
    ),
    LOverLobs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LOverLobs",
            coeffDict_,
            0.2
        )
    ),
    LOverMobs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LOverMobs",
            coeffDict_,
            0.0
        )
    ),
    bMin_(coeffDict_.get<scalar>("bMin")),
    flameFilteredG_(coeffDict_.get<bool>("flameFilteredG"))
 {}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

PDRkEpsilon::~PDRkEpsilon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool PDRkEpsilon::read()
{
    if (RASModel::read())
    {
        noTurbUntil_.readIfPresent(coeffDict_);
        LOverLobs_.readIfPresent(coeffDict_);
        LOverMobs_.readIfPresent(coeffDict_);
        return true;
    }

    return false;
}


void PDRkEpsilon::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        nut_ = Cmu_*sqr(k_)/epsilon_;
        nut_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    volScalarField divU("divU",fvc::div(phi_/fvc::interpolate(rho_)));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);

    volScalarField G(GName(), rho_*nut_*(tgradU() && dev(twoSymm(tgradU()))));

    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Add the blockage generation term so that it is included consistently
    // in both the k and epsilon equations
    const volScalarField& betav =
        U_.db().lookupObject<volScalarField>("betav");

    const volScalarField& Lobs =
        U_.db().lookupObject<volScalarField>("Lobs");

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    const volScalarField& b = mesh_.lookupObject<volScalarField>("b");

    const scalarField Cw(pow(mesh_.V(), 2.0/3.0));

    const PDRDragModel& drag =
        U_.db().lookupObject<PDRDragModel>("PDRDragModel");

    volScalarField GR(drag.Gk());

    volScalarField LD
    (
        "LD",
        (LOverLobs_)*(Lobs + dimensionedScalar("minLength", dimLength, VSMALL))
    );

    LD.primitiveFieldRef() = min(LD.primitiveField(), maxLOverCellW_*sqrt(Cw));

    const volScalarField LI(lCoef_*pow(k_, 3.0/2.0)/epsilon_);

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    ignition ign(combustionProperties, mesh_.time(), U.mesh());

    dimensionedVector ignLoc("ignLoc", dimLength, ign.sites()[0].location());

    dimensionedScalar filtRad2
    (
        "filtRad2",
        dimLength,
        6.0*ign.sites()[0].diameter()
    );

    if (flameFilteredG_)
    {
        const volScalarField filDist(mag(mesh_.C() - ignLoc));
        const volScalarField filterG(pos(filDist - filtRad2));
        const volScalarField filterdivU(pos(filDist - filtRad2));
        const volScalarField filterGR(pos(filDist - filtRad2)*pos(b - bMin_));

        G *= filterG;
        GR *= filterGR;
        divU *= filterdivU;
    }

    volScalarField Cl(C5_ + (C6_*((LI - LD)/LI)));
    Cl.max(0.0);

    if (mesh_.time() > noTurbUntil_)
    {
        tmp<fvScalarMatrix> epsEqn
        (
            betav*fvm::ddt(rho_, epsilon_)
          + fvm::div(phi_, epsilon_)
          - fvm::laplacian(rho_*DepsilonEff(), epsilon_)
        ==
            C1_*betav*G*epsilon_/k_
          + Cl*(epsilon_/k_)*GR
          - fvm::SuSp(((2.0/3.0)*C1_+C3_)*betav*rho_*divU, epsilon_)
          - fvm::Sp(C2_*betav*rho_*epsilon_/k_, epsilon_)
        );

        epsEqn.ref().relax();
        epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
        solve(epsEqn);
        bound(epsilon_, epsilonMin_);

        tmp<fvScalarMatrix> kEqn
        (
            betav*fvm::ddt(rho_, k_)
        + fvm::div(phi_, k_)
        - fvm::laplacian(rho_*DkEff(), k_)
        ==
          (betav*G + GR)
        - fvm::SuSp((2.0/3.0)*betav*rho_*divU, k_)
        - fvm::Sp(betav*rho_*epsilon_/k_, k_)
        );

        kEqn.ref().relax();
        solve(kEqn);
        bound(k_, kMin_);
    }
    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
