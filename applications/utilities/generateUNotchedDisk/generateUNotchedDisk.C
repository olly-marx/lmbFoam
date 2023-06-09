/*---------------------------------------------------------------------------*\
              Original work | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
              Modified work | Copyright (C) 2017-2018 Johan Roenby
-------------------------------------------------------------------------------

License
    This file is part of IsoAdvector, which is an unofficial extension to
    OpenFOAM.

    IsoAdvector is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IsoAdvector is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IsoAdvector. If not, see <http://www.gnu.org/licenses/>.

Application
    generateU

Description
    Generates velocity field for the classical test case with a notched disc
    in a rigid body rotation flow.

Author
    Johan Roenby, STROMNING, all rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const scalar pi = Foam::mathematicalConstant::pi;

    const scalarField x(mesh.C().component(vector::X));
    const scalarField y(mesh.C().component(vector::Y));
    const scalarField z(mesh.C().component(vector::Z));

    vectorField& Uc = U.internalField();
    Uc.replace(vector::X, -2*pi*(z - 0.5));
    Uc.replace(vector::Y, 0);
    Uc.replace(vector::Z, 2*pi*(x - 0.5));

    U.correctBoundaryConditions();

    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    // Calculating phi
    const vectorField& Cf(mesh.Cf().internalField());
    const scalarField Xf(Cf.component(vector::X));
    const scalarField Yf(Cf.component(vector::Y));
    const scalarField Zf(Cf.component(vector::Z));
    vectorField Uf(Xf.size());
    Uf.replace(0, -2*pi*(Zf - 0.5));
    Uf.replace(1, 0);
    Uf.replace(2, 2*pi*(Xf - 0.5));

    scalarField& phic = phi.internalField();
    const vectorField& Sfc = mesh.Sf().internalField();
    phic = Uf & Sfc;

    surfaceScalarField::GeometricBoundaryField& phibf = phi.boundaryField();
    const surfaceVectorField::GeometricBoundaryField& Sfbf =
        mesh.Sf().boundaryField();
    const surfaceVectorField::GeometricBoundaryField& Cfbf =
        mesh.Cf().boundaryField();

    forAll(phibf, patchi)
    {
        scalarField& phif = phibf[patchi];
        const vectorField& Sff = Sfbf[patchi];
        const vectorField& Cff = Cfbf[patchi];
        const scalarField xf(Cff.component(vector::X));
        const scalarField yf(Cff.component(vector::Y));
        const scalarField zf(Cff.component(vector::Z));
        vectorField Uf(xf.size());
        Uf.replace(0, -2*pi*(zf - 0.5));
        Uf.replace(1, 0);
        Uf.replace(2, 2*pi*(xf - 0.5));
        phif = Uf & Sff;
    }

    Info<< "Writing U and phi\n" << endl;
    phi.write();
    U.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
