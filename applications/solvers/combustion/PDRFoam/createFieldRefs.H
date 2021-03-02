const volScalarField& psi = thermo.psi();
scalar pMin = thermo.getOrDefault<scalar>("pMin", 0);
