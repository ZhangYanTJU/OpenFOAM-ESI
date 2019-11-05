
fvScalarMatrix paEqn
(
    fvm::d2dt2(pa) - sqr(c0)*fvc::laplacian(pa)
);

paEqn.relax();
paEqn.solve();

//pAcoustic.boundaryFieldRef().updateCoeffs();
