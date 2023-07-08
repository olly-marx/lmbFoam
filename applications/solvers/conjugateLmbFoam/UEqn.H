    // Note: using "smeared" effective viscosity - proper tangential stress
    // balance is not taken into account in this implementation of the GFM
    surfaceScalarField nuEff
    (
        "nuEff",
        twoPhaseProperties.nuf()
      + fvc::interpolate(turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nuEff, U)
      - (fvc::grad(U) & fvc::grad(nuEff))
     ==
        alpha1*bodyForce
    );

    UEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor
        (
            U.select(pimple.finalIter())
        )
    );

    solve
    (
        UEqn
      + fvc::reconstruct
        (
            // snGrad handles 1/rho jump from intFvDataVOF
            fvc::snGrad(pd)*mesh.magSf()
        ),
        mesh.solutionDict().solver(U.select(pimple.finalIter()))
    );
