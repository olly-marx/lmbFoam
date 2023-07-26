    // Now must define the boundary conditions for the static vector potential A
    const scalar mu0over4pi = 1.0e-7;

    Info<< "\nNow parallel\n" << endl;

    // This is where you define the global fields that you can access from all processors
    vectorField rGlobal = mesh.C();
    List<Field<vector>> allValuesR(Pstream::nProcs());
    allValuesR[Pstream::myProcNo()] = rGlobal;
    Pstream::gatherList(allValuesR);
    Pstream::scatterList(allValuesR);
    rGlobal =
    ListListOps::combine<Field<vector> >
    (
        allValuesR,
        accessOp<Field<vector> >()
    );

    scalarField vGlobal = mesh.V();
    List<Field<scalar>> allValuesV(Pstream::nProcs());
    allValuesV[Pstream::myProcNo()] = vGlobal;
    Pstream::gatherList(allValuesV);
    Pstream::scatterList(allValuesV);
    vGlobal =
    ListListOps::combine<Field<scalar> >
    (
        allValuesV,
        accessOp<Field<scalar> >()
    );

    vectorField J0Global = J0.internalField();
    List<Field<vector>> allValuesJ0(Pstream::nProcs());
    allValuesJ0[Pstream::myProcNo()] = J0Global;
    Pstream::gatherList(allValuesJ0);
    Pstream::scatterList(allValuesJ0);
    J0Global =
    ListListOps::combine<Field<vector> >
    (
        allValuesJ0,
        accessOp<Field<vector> >()
    );

    vectorField jGlobal = j.internalField();
    List<Field<vector>> allValuesj(Pstream::nProcs());
    allValuesj[Pstream::myProcNo()] = jGlobal;
    Pstream::gatherList(allValuesj);
    Pstream::scatterList(allValuesj);
    jGlobal =
    ListListOps::combine<Field<vector> >
    (
        allValuesj,
        accessOp<Field<vector> >()
    );



    // DEBUGGING

    //~ // Test rGlobal
    //~ vector testSum = vector::zero;
    //~ forAll(rGlobal, cellI)
    //~ {
        //~ testSum += rGlobal[cellI];
    //~ }
    //~ Info << "Sum over rGlobal = " << testSum << endl;

    //~ testSum = vector::zero;
    //~ forAll(mesh.C(), cellI)
    //~ {
        //~ testSum += mesh.C()[cellI];
    //~ }
    //~ Info << "Sum over rPar = " << testSum << endl;

    //~ // Test vGlobal
    //~ scalar testSumS = 0.0;
    //~ forAll(vGlobal, cellI)
    //~ {
        //~ testSumS += vGlobal[cellI];
    //~ }
    //~ Info << "Sum over vGlobal = " << testSumS << endl;

    //~ testSumS = 0.0;
    //~ forAll(mesh.C(), cellI)
    //~ {
        //~ testSumS += mesh.V()[cellI];
    //~ }
    //~ Info << "Sum over vPar = " << testSumS << endl;

    //~ // Test jGlobal
    //~ testSum = vector::zero;
    //~ forAll(jGlobal, cellI)
    //~ {
        //~ testSum += jGlobal[cellI];
    //~ }
    //~ Info << "Sum over jGlobal = " << testSum << endl;

    //~ testSum = vector::zero;
    //~ forAll(J0, cellI)
    //~ {
        //~ testSum += J0.internalField()[cellI];
    //~ }
    //~ Info << "Sum over jPar = " << testSum << endl;

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
	        
                forAll(rGlobal, cellJ)
                {
                    const vector& rJ = rGlobal[cellJ];
                    const vector r = rFace - rJ;
                    const scalar magR = mag(r);

                    A0Face += mu0over4pi*(J0Global[cellJ] / magR) * vGlobal[cellJ];
                    aFace += mu0over4pi*(jGlobal[cellJ] / magR) * vGlobal[cellJ];
                }

	        A0.boundaryField()[patchI][faceI] = A0Face;
	        a.boundaryField()[patchI][faceI] = aFace;
            }
	}
    }
