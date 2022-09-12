Info<< "\nConstructing " << CLOUD_BASE_TYPE_NAME << " cloud" << endl;
reactingTypeCloud parcels
(
    word(CLOUD_BASE_TYPE_NAME) & "Cloud1",
    g,
    rho,
    U,
    slgThermo
);
