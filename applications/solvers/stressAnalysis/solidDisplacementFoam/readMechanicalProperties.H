    Info<< "Reading mechanical properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dictionary& rhoDict(mechanicalProperties.subDict("rho"));
    word rhoType(rhoDict.get<word>("type"));

    autoPtr<volScalarField> rhoPtr;

    IOobject rhoIO
    (
        "rho",
        Time::timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (rhoType == "uniform")
    {
        scalar rhoValue(rhoDict.get<scalar>("value"));

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh,
                dimensionedScalar
                (
                    "rho",
                    dimMass/dimVolume,
                    rhoValue
                )
            )
        );
    }
    else if (rhoType == "field")
    {
        rhoIO.readOpt(IOobject::MUST_READ);

        rhoPtr.reset
        (
            new volScalarField
            (
                rhoIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for rho"
            << abort(FatalError);
    }

    volScalarField& rho = rhoPtr();

    const dictionary& EDict(mechanicalProperties.subDict("E"));
    word EType(EDict.get<word>("type"));

    autoPtr<volScalarField> EPtr;

    IOobject EHeader
    (
        "E",
        Time::timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (EType == "uniform")
    {
        scalar rhoEValue(EDict.get<scalar>("value"));

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh,
                dimensionedScalar
                (
                    "Erho",
                    dimMass/dimLength/sqr(dimTime),
                    rhoEValue
                )
            )
        );
    }
    else if (EType == "field")
    {
        EHeader.readOpt(IOobject::MUST_READ);

        EPtr.reset
        (
            new volScalarField
            (
                EHeader,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for E"
            << abort(FatalError);
    }

    volScalarField& rhoE = EPtr();

    autoPtr<volScalarField> nuPtr;

    IOobject nuIO
    (
        "nu",
        Time::timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& nuDict(mechanicalProperties.subDict("nu"));
    word nuType(nuDict.get<word>("type"));

    if (nuType == "uniform")
    {
        scalar nuValue(nuDict.get<scalar>("value"));
        nuPtr.reset
        (
            new volScalarField
            (
                nuIO,
                mesh,
                dimensionedScalar
                (
                    "nu",
                    dimless,
                    nuValue
                )
            )
        );
    }
    else if (nuType == "field")
    {
        nuIO.readOpt(IOobject::MUST_READ);
        nuPtr.reset
        (
            new volScalarField
            (
                nuIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for nu"
            << abort(FatalError);
    }

    volScalarField& nu = nuPtr();

    Info<< "Normalising E : E/rho\n" << endl;
    volScalarField E(rhoE/rho);

    Info<< "Calculating Lame's coefficients\n" << endl;

    volScalarField mu(E/(2.0*(1.0 + nu)));
    volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
    volScalarField threeK(E/(1.0 - 2.0*nu));

    if (mechanicalProperties.get<bool>("planeStress"))
    {
        Info<< "Plane Stress\n" << endl;

        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
        threeK = E/(1.0 - nu);
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }
