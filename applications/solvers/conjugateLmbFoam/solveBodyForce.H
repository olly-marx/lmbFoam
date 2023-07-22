{
    // Solve electric field equation for the constant applied voltage
    conductivity = alpha1*conductivity1 + (scalar(1) - alpha1)*conductivity2;
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

    // Calculate current density from the constant applied voltage
    J0 = 
	    conductivity*
	    (
	    - fvc::grad(Voltage)
	    );
    J0.correctBoundaryConditions();

    solidJ0 = 
	    solidConductivity*
	    (
	    - fvc::grad(solidVoltage)
	    );
    solidJ0.correctBoundaryConditions();

    // Solve electric field equation for the induced voltage in the fluid

    //while (simpleSolid.correctNonOrthogonal())
    //{
    //    coupledFvScalarMatrix inducedVoltageEqns(2);

    //    volScalarField DivUxB = (U & fvc::curl(Bext)) - (Bext & fvc::curl(U));

    //    // Fluid side
    //    fvScalarMatrix* inducedVoltageEqn =
    //        new fvScalarMatrix
    //        (
    //    	fvm::laplacian(conductivity,inducedVoltage)
    //    	== conductivity*DivUxB
    //        );

    //    inducedVoltageEqn->relax();

    //    // Solid side
    //    fvScalarMatrix* solidInducedVoltageEqn =
    //        new fvScalarMatrix
    //        (
    //            fvm::laplacian(solidConductivity, solidInducedVoltage)
    //        );

    //    solidInducedVoltageEqn->relax();

    //    // Add fluid equation
    //    inducedVoltageEqns.set(0, inducedVoltageEqn);

    //    // Add solid equation
    //    inducedVoltageEqns.set(1, solidInducedVoltageEqn);

    //    //inducedVoltageEqns.solve();
    //}

    volScalarField DivUxB = (U & fvc::curl(Bext)) - (Bext & fvc::curl(U));
    
    // Fluid side
    solve(
        fvm::laplacian(conductivity,inducedVoltage)
        == conductivity*DivUxB
    );

    // Calculate induced current density
    j = 
	    conductivity*
	    (
	    - fvc::grad(inducedVoltage)
	    + (U ^ Bext)
	    );
    j.correctBoundaryConditions();

    // The overall current density is the sum of the two
    J = ( J0 
	+ j
	);
    
    // Electric current field
    bodyForce = 
	    (J ^ Bext)
	    /(alpha1*rho1 + (1 - alpha1)*rho2);
    bodyForce.correctBoundaryConditions();
}
