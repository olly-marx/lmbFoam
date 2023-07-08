{
    // STEP 1: set up flux and velocity for Crank-Nicolson scheme

    // Initialize placeholders for current U and phi fields before the solution.
    // Note: needed for isoAdvector that uses truncated Crank-Nicolson scheme
    // (only flux and velocity fields are discretised with C-N scheme, while
    // alpha is always explicitly discretised).
    autoPtr<surfaceScalarField> curPhiPtr;
    autoPtr<volVectorField> curUPtr;

    // If topo changes have occured, we must not take the old flux
    if (!meshChanged)
    {
        Info<< "Setting up velocity and flux fields for isoAdvector." << endl;

        // Store current fields
        curPhiPtr.set(new surfaceScalarField("tmpPhi", phi));
        curUPtr.set(new volVectorField("tmpU", U));

        // Get Crank-Nicolson off-centre coefficient from isoAdvector subdict
        // Note: different definition as in ddtSchemes
        // - 0 = Euler explicit
        // - 0.5 = Full Crank-Nicolson
        // - 1 = Euler implicit
        const scalar CNCoeff = readScalar
        (
            mesh.solutionDict().subDict("isoAdvector").lookup
            (
                "CrankNicolsonCoeff"
            )
        );

        // Set current fields for Crank-Nicolson scheme
        phi = CNCoeff*phi + (1.0 - CNCoeff)*phi.oldTime();
        U = CNCoeff*U + (1.0 - CNCoeff)*U.oldTime();
    }
    else
    {
        Info<< "Using Euler implicit scheme for isoAdvector since we have"
            << " topological changes..." << endl;
    }


    // STEP 2: Advect the free surface using isoAdvector

    // Create isoAdvection object
    isoAdvection advector(alpha1, phi, U);

    // Advect the free surface
    advector.advect();


    // STEP 3: Restore flux and velocity fields after isoAdvector
    if
    (
        curPhiPtr.valid()
     && curUPtr.valid()
     && !meshChanged
    )
    {
        Info<< "Restoring velocity and flux fields after isoAdvector." << endl;

        phi = curPhiPtr();
        U = curUPtr();

        // Clear out stored fields
        curPhiPtr.clear();
        curUPtr.clear();
    }
    else if
    (
        (!curPhiPtr.valid() || !curUPtr.valid())
     && !meshChanged
    )
    {
        // Sanity check
        FatalError<< "Current flux and velocity fields are invalid after the"
                  << " isoAdvector advection step."
                  << nl
                  << "Something went wrong."
                  << exit(FatalError);
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;

    interface.correct();
    twoPhaseProperties.correct();

    // Update interfaceFvDataVOF object after obtaining new alpha field
    intFvDataVOF.update();
}
