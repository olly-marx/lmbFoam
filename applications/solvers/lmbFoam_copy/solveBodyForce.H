{
    // Solve electric field equation
    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;
    conductivity.correctBoundaryConditions();
    
    // Solve the Voltage equation which is the laplacian of the Voltage == the
    // divergence of U x B
    solve(fvm::laplacian(conductivity, Voltage)
		    //== fvc::div(U ^ Bext)
    );

    // Current density
    // OJM40
    currentDensity =
        conductivity *
	(
	- fvc::grad(Voltage)
	//+ (U ^ Bext)
	);
    currentDensity.correctBoundaryConditions();

    // Electric current field
    bodyForce = 
	(currentDensity ^ Bext)
	/ (alpha1*rho1 + (1 - alpha1)*rho2);	    // 1/rho
    bodyForce.correctBoundaryConditions();
}
