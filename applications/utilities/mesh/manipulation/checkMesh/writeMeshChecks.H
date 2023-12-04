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
    { writeChecksFormatType::JSON, "JSON" },
};

writeChecksFormatType writeChecksFormat(writeChecksFormatType::none);

auto writeMeshChecks = [](const fvMesh& mesh, const writeChecksFormatType fmt)
{
    if (Pstream::master())
    {
        switch (fmt)
        {
            case writeChecksFormatType::dictionary:
            {
                OFstream os(mesh.time().globalPath()/"checkMesh.dict");

                IOdictionary data
                (
                    IOobject
                    (
                        mesh.time().globalPath()/"checkMesh.dict",
                        mesh,
                        IOobject::NO_READ
                    )
                );

                data.writeHeader(os);

                mesh.data().meshDict().write(os, false);

                IOobject::writeEndDivider(os);

                Info<< "Writing mesh data to " << data.objectPath()
                    << nl << endl;

                break;
            }
            case writeChecksFormatType::JSON:
            {
                OFstream os(mesh.time().globalPath()/"checkMesh.json");

                Info<< "Writing mesh data to " << os.name() << nl << endl;

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
