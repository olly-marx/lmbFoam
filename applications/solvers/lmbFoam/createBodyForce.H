    Info<< "Reading field Voltage\n" << endl;
    volScalarField Voltage
    (
        IOobject
        (
            "Voltage",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Read external magnetic field
    volVectorField Bext
    (
        IOobject
        (
            "Bext",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Read induced magnetic field
    volVectorField B
    (
	IOobject
	(
	    "B",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );

    dimensionedScalar conductivity1
    (
        twoPhaseProperties.subDict("phase1").lookup("conductivity")
    );
    dimensionedScalar conductivity2
    (
        twoPhaseProperties.subDict("phase2").lookup("conductivity")
    );

    // Create electric conductivity field
    volScalarField conductivity
    (
        IOobject
        (
            "conductivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2
    );

    volVectorField bodyForce
    (
        IOobject
        (
            "bodyForce",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (conductivity*fvc::grad(Voltage) ^ Bext)/
        (alpha1*rho1 + (1 - alpha1)*rho2)
    );

    volVectorField currentDensity
    (
        IOobject
        (
            "currentDensity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2)*fvc::grad(Voltage)
    );













