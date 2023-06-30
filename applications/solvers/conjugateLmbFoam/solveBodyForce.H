{
    // Solve electric field equation
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
                fvm::laplacian(conductivity, Voltage)
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
        

    // Currect density
    currentDensity = conductivity*fvc::grad(Voltage);
    currentDensity.correctBoundaryConditions();

    solidCurrentDensity = solidConductivity*fvc::grad(solidVoltage);
    solidCurrentDensity.correctBoundaryConditions();
    
    // Electric current field
    bodyForce = (currentDensity ^ Bext)/(alpha1*rho1 + (1 - alpha1)*rho2);
    bodyForce.correctBoundaryConditions();
}
