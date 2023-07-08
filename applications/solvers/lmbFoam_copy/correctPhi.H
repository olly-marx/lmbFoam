    // Calculate the pressure correction equation with jump conditions, updating
    // the pressure as well. Necessary for topological changes

    // Update interface data since the topology is changing and we have a new
    // set of interface faces. Note: need to update old interface data as well
    Info<< "Updating interface data for topology change..." << endl;
    intFvDataVOF.update();
    intFvDataVOF.updateOld();

    // Note: reset dynamic pressure to zero in order to avoid using
    // non-orthogonal correction since the explicit pressure gradient depends
    // on old interface data (before the topo change)
    pd = dimensionedScalar("zero", pd.dimensions(), 0);

    // Calculate new pressure adhering to the new topology
#   include "pdInitEqn.H"

