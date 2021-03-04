    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiuReactionThermo> pThermo
    (
        psiuReactionThermo::New(mesh)
    );
    psiuReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "ha", "ea");
    basicMultiComponentMixture& composition = thermo.composition();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    volScalarField rhou
    (
        IOobject
        (
            "rhou",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rhou()
    );
    volScalarField muu
    (
        IOobject
        (
            "muu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.muu()
    );

    volScalarField& p = thermo.p();

    volScalarField& b = composition.Y("b");

    volScalarField combustFlag
    (
        IOobject
        (
            "combustFlag",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading fields betav, betai\n" << endl;
    volScalarField betav
    (
        IOobject
        (
            "betav",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField betai
    (
        IOobject
        (
            "betai",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Bv\n" << endl;
    volSymmTensorField Bv
    (
        IOobject
        (
            "Bv",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading fields Nv,nsv\n" << endl;
    volScalarField Nv
    (
        IOobject
        (
            "Nv",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField nsv
    (
        IOobject
        (
            "nsv",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading fields Blong\n" << endl;
    volSymmTensorField Blong
    (
        IOobject
        (
            "Blong",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "\nCreating field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"
    mesh.setFluxRequired(p.name());

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::New<compressible::RASModel>
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    #include "createDpdt.H"

    #include "createK.H"

    Info<< "Creating the unstrained laminar flame speed\n" << endl;
    autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed
    (
        laminarFlameSpeed::New(thermo)
    );

    Info<< "Reading strained laminar flame speed field Su\n" << endl;
    volScalarField Su
    (
        IOobject
        (
            "Su",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field Lobs\n" << endl;
    volScalarField Lobs
    (
        IOobject
        (
            "Lobs",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field CT\n" << endl;
    volSymmTensorField CT
    (
        IOobject
        (
            "CT",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    IOdictionary PDRProperties
    (
        IOobject
        (
            "PDRProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    //- Create the drag model
    autoPtr<PDRDragModel> drag = PDRDragModel::New
    (
        PDRProperties,
        turbulence,
        rho,
        U,
        phi
    );

    //- Create the flame-wrinkling model
    autoPtr<XiModel> flameWrinkling = XiModel::New
    (
        PDRProperties,
        thermo,
        turbulence,
        Su,
        rho,
        b,
        phi
    );

    Info<< "Calculating turbulent flame speed field St\n" << endl;
    volScalarField St
    (
        IOobject
        (
            "St",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        flameWrinkling->Xi()*Su
    );

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    if (composition.contains("ft"))
    {
        fields.add(composition.Y("ft"));
    }

    fields.add(b);
    fields.add(thermo.he());
    fields.add(thermo.heu());
    flameWrinkling->addXi(fields);

    #include "createMRF.H"
    #include "createFvOptions.H"

    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - b
    );
