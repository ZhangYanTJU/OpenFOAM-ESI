tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,ft_b_h_hu)")
    )
);

rhou = thermo.rhou();
muu = thermo.muu();

volScalarField Db("Db", turbulence->mut());
const volScalarField DbThermal
(
    "DbAlpha",
    turbulence->alphaEff()
    //thermo.alpha() + turbulence->mut()/Prt
);

if (ign.ignited())
{
    // Calculate the unstrained laminar flame speed
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Su = unstrainedLaminarFlameSpeed()();

    // Wrinkling due to turbulence
    const volScalarField& Xi = flameWrinkling->Xi();
    // Wrinkling due to obstacles
    const volScalarField& Xp = flameWrinkling->Xp();

    // progress variable
    // ~~~~~~~~~~~~~~~~~
    c = (1.0 - b);


    // Calculate flame normal etc.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    volVectorField n(fvc::reconstruct(fvc::snGrad(b)*mesh.magSf()));

    volScalarField mgb("mgb", mag(n));
    volScalarField ddtb("ddtb",fvc::ddt(b));

    dimensionedScalar dMgb("dMgb", mgb.dimensions(), SMALL);

    {
        volScalarField bc(b*c);

        dMgb += 1.0e-8*
            (bc*mgb)().weightedAverage(mesh.V())
           /(bc.weightedAverage(mesh.V()) + SMALL);
    }

    mgb += dMgb;

    surfaceVectorField Sfhat(mesh.Sf()/mesh.magSf());
    surfaceVectorField nfVec(fvc::interpolate(n));
    nfVec += Sfhat*(fvc::snGrad(b) - (Sfhat & nfVec));
    nfVec /= (mag(nfVec) + dMgb);
    surfaceScalarField nf("nf", mesh.Sf() & nfVec);
    n /= mgb;


    #include "StCorr.H"

    // Calculate turbulent flame speed flux
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    surfaceScalarField phiSt("phiSt", fvc::interpolate(rhou*StCorr*St)*nf);

    #include "StCourantNo.H"

    Db = flameWrinkling->Db();

    fvScalarMatrix bEqn
    (
        betav*fvm::ddt(rho, b)
      + mvConvection->fvmDiv(phi, b)
      + fvm::div(phiSt, b)
      - fvm::Sp(fvc::div(phiSt), b)
      - fvm::laplacian(Db, b)
     ==
        betav*fvOptions(rho, b)
        // + betav*rhou*St*StCorr*mag(fvc::grad(b))
    );


    // Add ignition cell contribution to b-equation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "ignite.H"

    bEqn.relax();
    fvOptions.constrain(bEqn);
    bEqn.solve();

    fvOptions.correct(b);

    Info<< "min(b) = " << min(b).value() << nl
        << "max(b) = " << max(b).value() << endl;

    if (composition.contains("ft"))
    {
        volScalarField& ft = composition.Y("ft");

        Info<< "Combustion progress = "
            << 100*(1.0 - b)().weightedAverage(mesh.V()*ft).value() << "%"
            << endl;
    }
    else
    {
        Info<< "Combustion progress = "
            << 100*(1.0 - b)().weightedAverage(mesh.V()).value() << "%"
            << endl;
    }


    // Correct the flame-wrinkling, smoothing
    flameWrinkling->correct(mvConvection);
    St = Xp*Xi*Su*combustFlag;

    if (PDRProperties.getOrDefault("smoothSt", false))
    {
        const scalar coeff = PDRProperties.get<scalar>("smoothStCoeff");
        fvc::smooth(St, coeff);
    }
}
