    // Relative momentum predictor

    tmp<fvVectorMatrix> tUrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevRhoReff(Urel)
      + rho*SRF->Su()
     ==
        fvOptions(rho, Urel)
    );
    fvVectorMatrix& UrelEqn = tUrelEqn.ref();

    UrelEqn.relax();

    fvOptions.constrain(UrelEqn);

    if (simple.momentumPredictor())
    {
        solve(UrelEqn == -fvc::grad(p));

        fvOptions.correct(Urel);
    }
