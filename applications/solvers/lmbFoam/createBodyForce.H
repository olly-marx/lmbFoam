    Info<< "Reading field Electric Potential\n" << endl;
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

    // Create electric conductivity field sigmaE
    volScalarField sigmaE
    (
        IOobject
        (
            "sigmaE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*sigmaE1 + (scalar(1) - alpha1)*sigmaE2
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
        sigmaE*
	(
	- fvc::grad(Voltage)
	+ (U ^ Bext)
	)
    );

    volVectorField fL
    (
        IOobject
        (
            "fL",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (J ^ Bext)/
        (alpha1*rho1 + (1 - alpha1)*rho2)
    );













