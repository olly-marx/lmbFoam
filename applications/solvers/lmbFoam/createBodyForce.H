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
	conductivity*
	(
	fvc::grad(Voltage)
	+ (U ^ Bext)
	)
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
        //(conductivity*fvc::grad(Voltage) ^ Bext)/
        (currentDensity ^ Bext) / // OJM40
        (alpha1*rho1 + (1 - alpha1)*rho2)
    );













