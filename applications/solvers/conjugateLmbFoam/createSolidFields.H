    Info<< "Reading field solidVoltage\n" << endl;
    volScalarField solidVoltage
    (
        IOobject
        (
            "solidVoltage",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    //Info<< "Reading field solidInducedVoltage\n" << endl;
    //volScalarField solidInducedVoltage
    //(
    //    IOobject
    //    (
    //        "solidInducedVoltage",
    //        runTime.timeName(),
    //        solidMesh,
    //        IOobject::MUST_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    solidMesh
    //);

    Info<< "Reading field solidConductivity\n" << endl;
    volScalarField solidConductivity
    (
        IOobject
        (
            "solidConductivity",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    // Note: mesh detached.  Do not evaluate current density
    volVectorField solidJ0
    (
        IOobject
        (
            "solidJ0",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, -2, 0, 0, 0, 1, 0),
            vector::zero
        )
    );
