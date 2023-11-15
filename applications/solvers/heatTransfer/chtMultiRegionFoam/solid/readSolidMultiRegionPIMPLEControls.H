const dictionary& pimple =
    mesh.solution().solutionDict("PIMPLE");

int nNonOrthCorr =
    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
