const dictionary& stressControl =
    mesh.solution().solutionDict("stressAnalysis");

bool compactNormalStress(stressControl.get<bool>("compactNormalStress"));
