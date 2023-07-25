    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
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

    // Electric potential field

#   include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "alpha1");

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    // OJM40

    dimensionedScalar nu1
    (
	twoPhaseProperties.subDict("phase1").lookup("nu")
    );
    dimensionedScalar nu2
    (
     	twoPhaseProperties.subDict("phase2").lookup("nu")
    );
    dimensionedScalar conductivity1
    (
     	twoPhaseProperties.subDict("phase1").lookup("conductivity")
    );
    dimensionedScalar conductivity2
    (
     	twoPhaseProperties.subDict("phase2").lookup("conductivity")
    );

    dimensionedScalar mu
    (
     	twoPhaseProperties.lookup("mu")
    );

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", pd.dimensions(), 0)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(pd, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    // Set required fluxes
    mesh.schemesDict().setFluxRequired(pd.name());
    mesh.schemesDict().setFluxRequired(alpha1.name());

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    // Create interfaceFvDataVOF object needed for Ghost Fluid Method
    interfaceFvDataVOF intFvDataVOF(mesh, alpha1);

    // Check whether the schemes are consistent with GFM. This solver uses
    // laplacian scheme for pressure equation, snGrad scheme for momentum
    // equation and gradient snGradCorr scheme for gradient in the
    // non-orthogonal correction
    intFvDataVOF.checkSchemes
    (
        "laplacian(rAU,pd)",
        "snGrad(pd)",
        "snGradCorr(pd)",
        true
    );

    // Initialise pd according to jump across the free surface
#   include "pdInitEqn.H"

    // Update total pressure field after initialising pd
#   include "pEqn.H"
