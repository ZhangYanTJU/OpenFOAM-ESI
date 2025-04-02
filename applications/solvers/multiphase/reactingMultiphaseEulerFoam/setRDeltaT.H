{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.2)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    surfaceScalarField maxPhi("maxPhi", phi);

    forAll(phases, phasei)
    {
        maxPhi = max(maxPhi, mag(phases[phasei].phi()));
    }

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(maxPhi)()()
       /((2*maxCo)*mesh.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    auto limits = gMinMax(rDeltaT.primitiveField());
    limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

    Info<< "Flow time scale min/max = "
        << limits.min() << ", " << limits.max() << endl;
}
