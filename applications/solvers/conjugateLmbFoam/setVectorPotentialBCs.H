    // Now must define the boundary conditions for the static vector potential A
    const scalar mu0over4pi = 1.0e-7;

//    reduce(J0, sumOp<vector>());
//    reduce(j, sumOp<vector>());

    Info << "Setting boundary conditions for static vector potential A0 and induced vector potential a" << endl;
    forAll(mesh.boundary(), patchI)
    {
	// Only if the patch is one of the following:
	// 1. liquidToSteelInterface
	// 2. liquidtoLithiumInterface
	// 3. topSurface
	
	if (mesh.boundary()[patchI].name() == "liquidToSteelInterface" ||
	    mesh.boundary()[patchI].name() == "liquidToLithiumInterface" ||
	    mesh.boundary()[patchI].name() == "topSurface")
	{
            const fvPatch& patch = mesh.boundary()[patchI];
            Info << "Visiting patch " << patch.name() << endl;
            forAll(patch, faceI)
            {
                const vector& rFace = patch.Cf()[faceI];

	        vector A0Face = vector::zero;
	        vector aFace = vector::zero;
	        
                forAll(mesh.cells(), cellJ)
                {

                    const vector& rJ = mesh.C()[cellJ];
                    const vector r = rFace - rJ;
                    const scalar magR = mag(r);

                    A0Face += mu0over4pi*(J0[cellJ] / magR) * mesh.V()[cellJ];
                    aFace += mu0over4pi*(j[cellJ] / magR) * mesh.V()[cellJ];
                }

	        A0.boundaryField()[patchI][faceI] = A0Face;
	        a.boundaryField()[patchI][faceI] = aFace;
            }
	}
    }
