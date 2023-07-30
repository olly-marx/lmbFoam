{

#   include "solveCurrentDensity.H"

#   include "solveMagneticField.H"

    volVectorField J0crossb
    (
        IOobject
        (
	    "J0crossb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	J0 ^ b
    );

    volVectorField J0crossBext
    (
        IOobject
        (
	    "J0crossBext",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	J0 ^ Bext
    );

    volVectorField J0crossB0
    (
        IOobject
        (
	    "J0crossB0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	J0 ^ B0
    );

    volVectorField jcrossb
    (
        IOobject
        (
	    "jcrossb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	j ^ b
    );

    volVectorField jcrossB0
    (
        IOobject
        (
	    "jcrossB0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	j ^ B0
    );

    volVectorField jcrossBext
    (
        IOobject
        (
	    "jcrossBext",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	j ^ Bext
    );
    
    // Electric current field
    bodyForce = 
	    (
	     J0crossb
	     + J0crossBext
	     + J0crossB0
	     + jcrossb
	     + jcrossB0
	     + jcrossBext
	    )
	    /(alpha1*rho1 + (1 - alpha1)*rho2);
    bodyForce.correctBoundaryConditions();
}
