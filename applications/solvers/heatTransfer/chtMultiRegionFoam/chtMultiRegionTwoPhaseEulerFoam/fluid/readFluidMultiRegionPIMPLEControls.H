    const dictionary& pimpleDict =
        mesh.solution().solutionDict("PIMPLE");

    Switch faceMomentum
    (
        pimpleDict.getOrDefault<Switch>("faceMomentum", false)
    );

    int nEnergyCorrectors
    (
        pimpleDict.getOrDefault<int>("nEnergyCorrectors", 1)
    );

    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);
