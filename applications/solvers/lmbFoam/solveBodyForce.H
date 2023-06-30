{
    // Solve electric field equation
    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;
    conductivity.correctBoundaryConditions();
    
    solve(fvm::laplacian(conductivity, Voltage));

    // Currect density
    // OJM40
    currentDensity =
        conductivity*
	(
	fvc::grad(Voltage)
	+ (U ^ Bext)
	);
    currentDensity.correctBoundaryConditions();

    // Electric current field
    //bodyForce = (currentDensity ^ Bext)/(alpha1*rho1 + (1 - alpha1)*rho2);
    bodyForce = (currentDensity ^ Bext)/(alpha1*rho1 + (1 - alpha1)*rho2); // OJM40
    bodyForce.correctBoundaryConditions();
}
