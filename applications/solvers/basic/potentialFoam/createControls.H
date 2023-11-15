const dictionary& potentialFlow
(
    mesh.solution().solutionDict("potentialFlow")
);

const int nNonOrthCorr
(
    potentialFlow.getOrDefault<int>("nNonOrthogonalCorrectors", 0)
);
