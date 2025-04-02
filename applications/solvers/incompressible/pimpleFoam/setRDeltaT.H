{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.8)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi))()()
       /((2*maxCo)*mesh.V())
    );

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    {
        auto limits = gMinMax(rDeltaT.primitiveField());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

        Info<< "Flow time scale min/max = "
            << limits.min() << ", " << limits.max() << endl;
    }

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    {
        auto limits = gMinMax(rDeltaT.primitiveField());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

        Info<< "Smoothed flow time scale min/max = "
            << limits.min() << ", " << limits.max() << endl;
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        auto limits = gMinMax(rDeltaT.primitiveField());
        limits.reset(1/(limits.max()+VSMALL), 1/(limits.min()+VSMALL));

        Info<< "Damped flow time scale min/max = "
            << limits.min() << ", " << limits.max() << endl;
    }
}
