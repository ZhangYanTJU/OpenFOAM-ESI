fvSolution solutionDict(runTime);

const bool coupled(solutionDict.getOrDefault("coupledEnergyField", false));

autoPtr<lduPrimitiveMeshAssembly> meshPtr;
autoPtr<fvMatrixAssembly> fvMatrixAssemblyPtr;

if (coupled)
{
    UPtrList<fvMesh> meshes(fluidRegions.size() + solidRegions.size());

    label meshi = 0;

    for (fvMesh& mesh : fluidRegions)
    {
        meshes.set(meshi, &mesh);
        ++meshi;
    }

    for (fvMesh& mesh : solidRegions)
    {
        meshes.set(meshi, &mesh);
        ++meshi;
    }

    Info<< "Creating coupled mesh." << nl << endl;

    IOobject io
    (
        "assemblyLdu",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    );

    meshPtr.reset
    (
        new lduPrimitiveMeshAssembly(meshes, io)
    );

    Info<< "Creating coupled matrix." << nl << endl;

    fvMatrixAssemblyPtr.reset
    (
        new fvMatrixAssembly
        (
            meshPtr(),
            dimEnergy/dimTime,
            "he"
        )
    );

    Info<< "Created coupled matrix." << nl << endl;
}
