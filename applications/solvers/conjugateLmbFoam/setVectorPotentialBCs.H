    // Now must define the boundary conditions for the static vector potential A
    const scalar mu0over4pi = 1.0e-7;

    Info << "Setting boundary conditions for static vector potential A0 and induced vector potential a" << endl;
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        Info << "Visiting patch " << patch.name() << endl;
        forAll(patch, faceI)
        {
            const vector& rFace = patch.Cf()[faceI];

            forAll(mesh.cells(), cellJ)
            {

                const vector& rJ = mesh.C()[cellJ];
                const vector r = rFace - rJ;
                const scalar magR = mag(r);

                // Might need to modify which J you are using here

                A0.boundaryField()[patchI][faceI] += mu0over4pi*(J0[cellJ] / magR) * mesh.V()[cellJ];
                a.boundaryField()[patchI][faceI] += mu0over4pi*(j[cellJ] / magR) * mesh.V()[cellJ];
            }
        }
    }
