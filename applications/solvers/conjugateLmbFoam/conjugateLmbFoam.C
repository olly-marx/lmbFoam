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
    conjugateLmbFoam

Description
    Transient solver for 2 incompressible fluids with Ghost Fluid Method for
    treatment of jump conditions at the interface and geometric isoAdvector
    advection method based on Volume--of--Fluid approach.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    Dynamic mesh motion is executed within the PIMPLE loop. Intended to be used
    with adaptive mesh refinement.

    Added electric field model for liquid metal battery.
    Using conjugate fluid-solid solution for te electric field calculation.

Author
    Hrvoje Jasak, Wikki Ltd. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"
#include "regionCouplePolyPatch.H"
#include "interfaceFvDataVOF.H"
#include "isoAdvection.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "simpleControl.H"
#include "dynamicFvMesh.H"
#include "numericsFunctions.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createSolidMesh.H"

    pimpleControl pimple(mesh);
    simpleControl simpleSolid(solidMesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createControls.H"
#   include "createFields.H"
#   include "createSolidFields.H"
#   include "createMHDFields.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "alphaCourantNo.H"
#       include "setAlphaDeltaT.H"
#       include "setVectorPotentialBCs.H"

        // Detach coupled patches
#       include "detachPatches.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Mesh motion, free surface, and pressure-velocity corrector
        while (pimple.loop())
        {
            // Attached coupled patches
#           include "attachPatches.H"

            // Solve current to determine the electromagnetic force
#           include "solveBodyForce.H"

            // Detach coupled patches
#           include "detachPatches.H"

            // Correct mesh motion
#           include "correctMeshMotion.H"

            // Momentum predictor
#           include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pdEqn.H"
            }

            // Update total pressure field
#           include "pEqn.H"

            // Update free surface
#           include "alphaEqn.H"

	    // Check the lithium interface boundary patch for contact with the
	    // lower electrode, this occurs when the alpha field is greater than
	    // 0.5. In this case end the simulation.
	    //label patchi = mesh.boundaryMesh().findPatchID("lithiumInterface");
	    //if(patchi >= 0)
	    //{
	    //        scalarField alphaLithiumInterface = alpha1.boundaryField()[patchi];

	    //        scalar alphaMax = gMax(alphaLithiumInterface);
	    //        Info<< "alphaMax = " << alphaMax << endl;

	    //        if(alphaMax > 0.4)
	    //        {
	    //    	Info<< "FATAL CONDITION: SHORT CIRCUIT DETECTED" << endl;

	    //    	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    //    	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    //    	    << nl << "End\n" << endl;
	    //    	// Write data and then exit
	    //    	Info<< "Writing data at time " << runTime.timeName() << endl;
	    //    	alpha1.write();
	    //    	U.write();
	    //    	bodyForce.write();
	    //    	J.write();
	    //    	exit(0);
	    //        }
	    //        else if(alphaMax > 1.0e-10)
	    //        {
	    //    	Info<< "Nearing Short Circuit Condition" << nl 
	    //    		<<"Writing data at time " << runTime.timeName() << endl;
	    //    	alpha1.write();
	    //        }
	    //}

            turbulence->correct();

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
