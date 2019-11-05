
Info << "\nReading pa" << endl;

volScalarField pa
( 
	IOobject 
	( 
		"pa", 
        runTime.timeName(),
        mesh,
	    IOobject::MUST_READ,
        IOobject::AUTO_WRITE
	), 
	mesh
);
/*
Info << "\nCreating volume fields" << endl;

volScalarField w
(
    IOobject
    (
        "w",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("w",dimensionSet(0,1,0,0,0,0,0),0)
);

volScalarField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("v",dimensionSet(0,1,-1,0,0,0,0),0)
);

volScalarField a
(
    IOobject
    (
        "a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("a",dimensionSet(0,1,-2,0,0,0,0),0)
);

volScalarField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("F",dimensionSet(1,-1,-2,0,0,0,0),0)
);
*/
