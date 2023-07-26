{

#   include "solveCurrentDensity.H"

#   include "solveMagneticField.H"
    
    // Electric current field
    bodyForce = 
	    (J ^ B)
	    /(alpha1*rho1 + (1 - alpha1)*rho2);
    bodyForce.correctBoundaryConditions();
}
