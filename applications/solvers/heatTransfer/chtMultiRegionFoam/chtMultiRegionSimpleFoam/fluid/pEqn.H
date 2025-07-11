{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    //tUEqn.clear();

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);

    // Solve pressure
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p = p_rgh + rho*gh;

    #include "continuityErrs.H"

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume)
    {
        if (!compressible)
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
        }
        else
        {
            p += (initialMass - fvc::domainIntegrate(psi*p))
                /compressibility;
        }
        p_rgh = p - rho*gh;
    }

    rho = thermo.rho();

    rho.clamp_range(rhoMin[i], rhoMax[i]);
    rho.relax();

    Info<< "Min/max rho:" << min(rho).value() << ' '
        << max(rho).value() << endl;
}
