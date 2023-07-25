    //Info << "Setting B0 values" << endl;
    //forAll(B0, cellI)
    //{
    //    const vector& rI = mesh.C()[cellI];

    //    forAll(mesh.cells(), cellJ)
    //    {
    //        if (cellJ == cellI)
    //        {
    //            continue;
    //        }

    //        const vector& rJ = mesh.C()[cellJ];
    //        const vector r = rI - rJ;
    //        const scalar magR = mag(r);

    //        // Calculate the integrand
    //        const scalar mu0over4pi = 1.0e-7;
    //        B0[cellI] += (mu0over4pi*(J0[cellJ]^r)*mesh.V()[cellJ])/pow(magR,3);
    //    }
    //}

    // Now solve for the induced vector potential a
    volScalarField Dsigmamu = 1/(conductivity*mu);

    // Now solve for the static vector potential A by the induction equation
    solve(
	     fvm::laplacian(Dsigmamu,A0)
	     - fvc::grad(Voltage)
    );

    // Calculate the static magnetic field
    B0 = fvc::curl(A0);

    //*****************************************************************************************
    // Now solve for the induced magnetic field

    solve(
	     fvm::laplacian(Dsigmamu,a)
	     + (U ^ B)
	     + (U ^ (fvc::curl(a)))
	     - fvc::grad(inducedVoltage)
    );

    // Calculate the induced magnetic field
    b = fvc::curl(a);

    // Calculate the total magnetic field
    B = Bext + B0 + b;
