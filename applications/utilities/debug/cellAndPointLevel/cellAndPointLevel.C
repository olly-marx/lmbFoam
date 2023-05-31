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
    cellAndPointLevel

Description
    Write cellLevel field as volScalarField and pointLevelField as
    pointScalarField.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "labelIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Read cell level
    labelIOField cellLevel
    (
        IOobject
        (
            "cellLevel",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read point level
    labelIOField pointLevel
    (
        IOobject
        (
            "pointLevel",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create cell level as volScalarField
    volScalarField cellLevelField
    (
        IOobject
        (
            "cellLevelField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    scalarField& clIn = cellLevelField.internalField();
    forAll (clIn, cellI)
    {
        clIn[cellI] = cellLevel[cellI];
    }

    Info<< "Writing cellLevelField..." << endl;
    cellLevelField.write();

    // Create point level as pointScalarField
    const pointMesh pointM(mesh);

    pointScalarField pointLevelField
    (
        IOobject
        (
            "pointLevelField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointM,
        dimensionedScalar("zero", dimless, 0.0)
    );
    scalarField& plIn = pointLevelField.internalField();
    forAll (plIn, pointI)
    {
        plIn[pointI] = pointLevel[pointI];
    }

    Info<< "Writing pointLevelField..." << endl;
    pointLevelField.write();

    Info<< "End\n" << endl;

    return 0;
}
