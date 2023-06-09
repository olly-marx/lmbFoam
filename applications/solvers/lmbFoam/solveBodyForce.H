{
    // Solve electric field equation
    sigmaE = alpha1*sigmaE1 + (scalar(1) - alpha1)*sigmaE2;
    sigmaE.correctBoundaryConditions();
    
    solve(fvm::laplacian(sigmaE, Voltage)
    );

    // Currect density
    J =
        sigmaE *
	(
	- fvc::grad(Voltage)
	+ (U ^ Bext)
	);
    J.correctBoundaryConditions();

    // Electric current field
    fL = 
	(J ^ Bext)
	/(alpha1*rho1 + (1 - alpha1)*rho2);
    fL.correctBoundaryConditions();
}
