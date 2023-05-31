#   include "createTimeControls.H"

    // Mesh motion controls
    Switch correctPhi
    (
        pimple.dict().lookupOrDefault<Switch>("correctPhi", true)
    );

    Switch checkMeshCourantNo
    (
        pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", true)
    );

    scalar startMotionTime
    (
        pimple.dict().lookupOrDefault<scalar>("startMotionTime", 0.0)
    );
