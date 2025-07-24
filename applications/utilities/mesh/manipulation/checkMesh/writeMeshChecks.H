enum class writeChecksFormatType
{
    none,
    dictionary,
    JSON
};
const Enum<writeChecksFormatType> writeChecksFormatTypeNames
{
    { writeChecksFormatType::none, "none" },
    { writeChecksFormatType::dictionary, "dictionary" },
    { writeChecksFormatType::JSON, "json" },
};

writeChecksFormatType writeChecksFormat(writeChecksFormatType::none);

auto writeMeshChecks = [&](const fvMesh& mesh, const writeChecksFormatType fmt)
{
    if (UPstream::master())
    {
        fileName path(mesh.time().globalPath()/"checkMesh");
        if (mesh.name() != polyMesh::defaultRegion)
        {
            path += "_" + mesh.name();
        }

        OFstream os(path.ext(writeChecksFormatTypeNames[fmt]));

        Info<< "Writing mesh data to " << os.name() << nl << endl;

        switch (fmt)
        {
            case writeChecksFormatType::dictionary:
            {
                IOdictionary data
                (
                    IOobject
                    (
                        os.name(),
                        mesh,
                        IOobject::NO_READ
                    )
                );

                data.writeHeader(os);

                mesh.data().meshDict().write(os, false);

                IOobject::writeEndDivider(os);

                break;
            }
            case writeChecksFormatType::JSON:
            {
                JSONformatter json(os);
                json.writeDict(mesh.data().meshDict());

                break;
            }
            default:
            {
                // Do nothing
            }
        }
    }
};
