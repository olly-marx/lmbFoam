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
    initWaveField

Description
    Initialise alpha field to a hard-coded shape

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sphereDistance.H"
#include "ImmersedCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("resetPatchAlpha", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Read fields

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
    scalarField& alphaIn = alpha1.internalField();

    const scalarField& VIn = mesh.V().field();

    // Construct a distance object for the ImmersedCell construction
    sphereDistance dist(vector(0, 0, 0), 1e-4);

    Info<< "Setting alpha... ";
    forAll (alphaIn, cellI)
    {
        // Use the ImmersedCell class insted of wetCell (IG 7/Jan/2019)
        ImmersedCell<sphereDistance> wc(cellI, mesh, dist);
        //wetCell wc(cellI, mesh, waveT);

        if (wc.isAllWet())
        {
            alphaIn[cellI] = 1;
        }
        else if (wc.isAllDry())
        {
            alphaIn[cellI] = 0;
        }
        else
        {
            alphaIn[cellI] = wc.wetVolume()/VIn[cellI];
        }
    }
    Info<< " done." << nl << endl;

    if (args.optionFound("resetPatchAlpha"))
    {
        Info<< "Resetting patch alpha" << endl;

        forAll (alpha1.boundaryField(), patchI)
        {
            alpha1.boundaryField()[patchI] ==
                alpha1.boundaryField()[patchI].patchInternalField();
        }
    }
    else
    {
        alpha1.correctBoundaryConditions();
    }

    Info<< "Writing alpha1" << endl;
    alpha1.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
