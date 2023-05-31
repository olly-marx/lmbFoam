/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    dummyLoadBalance

Description
    Perform a single dummy load balancing step and create two fields: alpha1 and
    procIndex to easily visualize the mapping.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "loadBalanceFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    loadBalanceFvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Create alpha1 field for checking
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Create processor index field and write
    volScalarField procIndex
    (
        IOobject
        (
            "procIndex",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    );
    scalarField& procIndexIn = procIndex.internalField();

    forAll (procIndexIn, cellI)
    {
        procIndexIn[cellI] = Pstream::myProcNo();
    }

    procIndex.correctBoundaryConditions();
    procIndex.write();

    // Increment time, perform load balancing and write the new mesh
    runTime++;

    Info<< "Performing dummy load balancing step..." << endl;
    mesh.update();

    Info<< "Writing load balanced mesh for time " << runTime.value() << endl;
    mesh.write();

    // Write processor index for the new configuration
    forAll (procIndexIn, cellI)
    {
       procIndexIn[cellI] = Pstream::myProcNo();
    }
    procIndex.correctBoundaryConditions();
    procIndex.write();

    // Write alpha field for the new configuration
    alpha1.correctBoundaryConditions();
    alpha1.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
