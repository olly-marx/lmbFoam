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
      //- fvc::div(phiB, 2.0*DBU*(Bext + B)) // OJM40
      - fvm::laplacian(nuEff, U)
      - (fvc::grad(U) & fvc::grad(nuEff))
      //+ fvc::grad(DBU*magSqr(Bext + B)) // OJM40
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
