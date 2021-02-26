if (ign.ignited())
{
    volScalarField& heau = thermo.heu();

    fvScalarMatrix heauEqn
    (
        betav*fvm::ddt(rho, heau) + mvConvection->fvmDiv(phi, heau)
      + (betav*fvc::ddt(rho, K) + fvc::div(phi, K))*rho/thermo.rhou()
      + (
            heau.name() == "eau"
          ? fvc::div
            (
                phi/fvc::interpolate(rho),
                p,
                "div(phiv,p)"
            )*rho/thermo.rhou()
          : -betav*dpdt*rho/thermo.rhou()
        )
      - fvm::laplacian(DbThermal, heau)
     ==
        betav*fvOptions(rho, heau)
    );

    fvOptions.constrain(heauEqn);

    heauEqn.solve();

    fvOptions.correct(heau);
}
