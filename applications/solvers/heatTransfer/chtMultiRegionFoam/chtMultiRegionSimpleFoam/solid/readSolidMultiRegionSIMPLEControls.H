const dictionary& simple =
    mesh.solution().solutionDict("SIMPLE");

const int nNonOrthCorr =
    simple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
