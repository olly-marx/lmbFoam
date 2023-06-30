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

    Info<< "Reading field conductivity\n" << endl;
    volScalarField conductivity
    (
        IOobject
        (
            "conductivity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;

    // Note: mesh detached.  Do not evaluate body force
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
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -2, 0, 0, 0, 0),
            vector::zero
        )
    );

    // Note: mesh detached.  Do not evaluate current density
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
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, -2, 0, 0, 0, 1, 0),
            vector::zero
        )
    );
