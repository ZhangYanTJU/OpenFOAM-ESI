Info<< "\nConstructing coal cloud" << endl;
reactingMultiphaseCloud coalParcels
(
    "coalCloud1",
    g,
    rho,
    U,
    slgThermo
);

Info<< "\nConstructing limestone cloud" << endl;
thermoCloud limestoneParcels
(
    "limestoneCloud1",
    g,
    rho,
    U,
    slgThermo
);
