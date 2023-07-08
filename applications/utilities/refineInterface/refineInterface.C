/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    refineInterface

Description
    Refines a polyhedral mesh close to the interface as defined by the
    interfaceFvDataVOF object. Number of refinement iterations is the same as
    maxRefinementLevel from the refinement dictionary.

Authors
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.
    (Idea by Inno Gatin, FMENA Zagreb)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicPolyRefinementFvMesh.H"
#include "interfaceFvDataVOF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Create mesh manually
    dynamicPolyRefinementFvMesh mesh
    (
        IOobject
        (
            dynamicFvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Read alpha field for refinement
    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Read velcity field in order to map it
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create interfaceFvDataVOF object needed for refinement strategy
    // "interfaceProximityRefinement"
#   include "readGravitationalAcceleration.H"
    interfaceFvDataVOF intFvDataVOF(mesh, alpha1);

    // Get number of refinement levels, which is equal to number of iterations
    // that we need to refine (at least)
    const label nIters = 10;  // Hacked. HJ, 28/Apr/2023
        // readLabel(mesh.refinementDict_.lookup("maxRefinementLevel"));

    Info<< "Performing: " << nIters << " refinement iterations." << nl << endl;

    // Loop through time because dynamicPolyRefinementFvMesh::update() depends
    // on it
    label i = 0;

    while (runTime.loop())
    {
        Info<< "Refinement iteration: " << i << endl;

        if (i++ == nIters)
        {
            break;
        }

        mesh.update();

        // Update interface data since the topology is changing and we have a
        // new set of interface faces.
        Info<< "Updating interface data for topology change..." << endl;
        intFvDataVOF.update();
        intFvDataVOF.updateOld();
    }

    // Write the final mesh
    mesh.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
