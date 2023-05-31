bool meshChanged = false;

if (runTime.value() > startMotionTime)
{
    // Make the fluxes absolute
    fvc::makeAbsolute(phi, U);

    meshChanged = mesh.update();

#   include "volContinuity.H"

    if (checkMeshCourantNo)
    {
#        include "meshCourantNo.H"
    }

    // Mesh motion update
    if (correctPhi && meshChanged)
    {
#        include "correctPhi.H"
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    if (meshChanged)
    {
        // Calculate and report Courant number
#       include "CourantNo.H"

        // Update gh and ghf
        gh = (g & mesh.C());
        ghf = (g & mesh.Cf());
    }

    // Correct velocity boundary conditions after mesh motion, with relative
    // fluxes
    U.correctBoundaryConditions();
}
