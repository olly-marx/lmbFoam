    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;

    // Solve electric field equation for the static applied voltage
    conductivity.correctBoundaryConditions();
    solidConductivity.correctBoundaryConditions();

    while (simpleSolid.correctNonOrthogonal())
    {
        coupledFvScalarMatrix VoltageEqns(2);

        // Fluid side
        fvScalarMatrix* VoltageEqn =
            new fvScalarMatrix
            (
                fvm::laplacian(conductivity,Voltage)
            );

        VoltageEqn->relax();

        // Solid side
        fvScalarMatrix* solidVoltageEqn =
            new fvScalarMatrix
            (
                fvm::laplacian(solidConductivity, solidVoltage)
            );

        solidVoltageEqn->relax();

        // Add fluid equation
        VoltageEqns.set(0, VoltageEqn);

        // Add solid equation
        VoltageEqns.set(1, solidVoltageEqn);

        VoltageEqns.solve();
    }
    
    // Calculate static current density from the constant applied voltage
    J0 = 
	conductivity*
	(- fvc::grad(Voltage));
    J0.correctBoundaryConditions();

    solidJ0 =
    solidConductivity*
    (
    - fvc::grad(solidVoltage)
    );
    solidJ0.correctBoundaryConditions();

    volScalarField DivUxB = (U & fvc::curl(B)) - (B & fvc::curl(U));
    
    // Fluid side
    // Print statement to check that the equation is being solved
    solve(
        fvm::laplacian(conductivity,inducedVoltage)
        == conductivity*DivUxB
    );
    inducedVoltage.correctBoundaryConditions();

    // Calculate induced current density
    j = 
	    conductivity*
	    (
	    - fvc::grad(inducedVoltage)
	    + (U ^ Bext)
	    );
    j.correctBoundaryConditions();

    // The overall current density is the sum of the two
    J = J0 + j;
