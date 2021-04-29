Foam::PtrList<Foam::fvMesh> meshes;

Foam::word regionName = Foam::fvMesh::defaultRegion;

if (args.readIfPresent("region", regionName))
{
    Foam::Info
        << "Create mesh " << regionName << " for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    meshes.setSize(1);
    meshes.set
    (
        0,
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            ),
            false
        )
    );
}
else if (args.found("allRegions"))
{
    wordList regionNames(regionProperties(runTime).names());
    meshes.setSize(regionNames.size());
    forAll(regionNames, i)
    {
        Foam::Info
            << "Create mesh " << regionNames[i] << " for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;

        meshes.set
        (
            i,
            new Foam::fvMesh
            (
                Foam::IOobject
                (
                    regionNames[i],
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                ),
                false
            )
        );
    }
}
else
{
    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    meshes.setSize(1);
    meshes.set
    (
        0,
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            ),
            false
        )
    );
}


for (auto& mesh : meshes)
{
    mesh.init(true);   // initialise all (lower levels and current)
}
