    Info << "Setting B0 values" << endl;
    forAll(B0, cellI)
    {
        const vector& rI = mesh.C()[cellI];

        forAll(mesh.cells(), cellJ)
        {
            if (cellJ == cellI)
            {
                continue;
            }

            const vector& rJ = mesh.C()[cellJ];
            const vector r = rI - rJ;
            const scalar magR = mag(r);

            // Calculate the integrand
            const scalar mu0over4pi = 1.0e-7;
            B0[cellI] += (mu0over4pi*(J0[cellJ]^r)*mesh.V()[cellJ])/pow(magR,3);
        }
    }
