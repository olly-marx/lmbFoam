// Initialize pressure with deltaT diffusion coefficient
surfaceScalarField phiU
(
    "phiU",
    (fvc::interpolate(U) & mesh.Sf())
);

phi = phiU;

// Adjust the fluxes
adjustPhi(phi, U, pd);

while (pimple.correctNonOrthogonal())
{
    // Pressure equation where laplacian handles 1/rho jump from
    // intFvDataVOF.
    fvScalarMatrix pdInitEqn
    (
        fvm::laplacian
        (
            runTime.deltaT(), pd, "laplacian(rAU,pd)"
        )
     ==
        fvc::div(phi)
    );

    pdInitEqn.setReference(pRefCell, pRefValue);

    pdInitEqn.solve(mesh.solutionDict().solver(pd.select("pdFinal")));

    phi -= pdInitEqn.flux();
}

// No velocity reconstruction as there's no ap. VV, 3/Apr/2015
