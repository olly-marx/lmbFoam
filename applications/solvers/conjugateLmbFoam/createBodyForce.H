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

    Info<< "Reading field inducedVoltage\n" << endl;
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
    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;

    // Note: mesh detached.  Do not evaluate current density
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
