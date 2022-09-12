const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
kinematicCloud parcels
(
    kinematicCloudName,
    g,
    rhoInf,
    U,
    muc
);

