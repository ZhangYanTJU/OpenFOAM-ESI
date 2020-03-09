volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}
else
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
}


MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

if (runTime.outputTime())
{
    volScalarField divPhiHbyA("divPhiHbyA", fvc::div(phiHbyA));
    divPhiHbyA.write();
    HbyA.write();
    rAU.write();
}

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();

        if (runTime.outputTime())
        {
            word iterN("divPhi" + name(pimple.corrPISO()));
            volScalarField divPhi(iterN, fvc::div(phi));
            divPhi.write();
            volScalarField divMeshFlux("divMeshFlux",
                fvc::surfaceIntegrate(pEqn.flux()));

            //volScalarField divMeshFlux("divMeshFlux", fvc::div(pEqn.flux()));
            divMeshFlux.write();
        }
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

if (runTime.outputTime())
{
    surfaceScalarField phiTmp(fvc::flux(U));
    volScalarField divPhi("divUTemp", fvc::div(phiTmp));
    divPhi.write();
}

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
