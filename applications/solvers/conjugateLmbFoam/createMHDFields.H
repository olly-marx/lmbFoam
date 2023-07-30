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

    // Initialisation of the static fields
    volVectorField J0
    (
        IOobject
        (
            "J0",
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

    // Initialisation of the static vector potential field
    volVectorField A0
    (
	IOobject
	(
	    "A0",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	),
	mesh
    );

    // Solve the initial magnetic field based on the constant vector potential
    volVectorField B0
    (
        IOobject
        (
            "B0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
	dimensionedVector
	(
	    "zero",
	    dimensionSet(1, 0, -2, 0, 0, -1, 0),
	    vector::zero
	)
    );

    // Initialisation of the induced voltage field
    volScalarField inducedVoltage
    (
	IOobject
	(
	    "inducedVoltage",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );

    // Initialisation of the induced current density field
    volVectorField j
    (
        IOobject
        (
            "j",
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

    // Initialisation of the induced vector potential field
    volVectorField a
    (
	IOobject
	(
	    "a",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	),
	mesh
    );

    // Initialisation of the induced magnetic field
    volVectorField b
    (
        IOobject
        (
            "b",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -2, 0, 0, -1, 0),
            vector::zero
        )
    );

    // The calculated fields are not constant, so they are not calculated here
    volVectorField J
    (
        IOobject
        (
            "J",
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

    // Calculated magnetic field B
    volVectorField B
    (
	IOobject
	(
	    "B",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	B0 + Bext
    );

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

