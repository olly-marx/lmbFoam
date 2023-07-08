/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

#include "dynamicLoadBalancedRefinementFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug
(
    dynamicLoadBalancedRefinementFvMesh,
    0
);

addToRunTimeSelectionTable
(
    dynamicFvMesh,
    dynamicLoadBalancedRefinementFvMesh,
    IOobject
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancedRefinementFvMesh::dynamicLoadBalancedRefinementFvMesh
(
    const IOobject& io
)
:
    // Note: create base class with the name of this class to be able to define
    // all the controls in dynamicLoadBalancedRefinementFvMeshCoeffs
    // subdictionary in dynamicMeshDict (see dynamicPolyRefinementFvMesh
    // constructor). VV, 31/12/2018.
    dynamicPolyRefinementFvMesh(io, typeName),
    loadBalance_(refinementDict_.lookup("loadBalance")),
    imbalanceThreshold_
    (
        readScalar(refinementDict_.lookup("imbalanceThreshold"))
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicLoadBalancedRefinementFvMesh::update()
{
    // Refine the mesh
    bool hasChanged = dynamicPolyRefinementFvMesh::update();

    // Load balance only if:
    // 1. Parallel run,
    // 2. Load balance specified,
    // 3. AMR has been performed
    if
    (
        Pstream::parRun()
     && loadBalance_
     && hasChanged
    )
    {
        // Check whether we really need to perform load balancing by looking at
        // number of cells in each processors
        labelList nCellsPerProc(Pstream::nProcs());

        // Fill in my part of the list
        nCellsPerProc[Pstream::myProcNo()] = this->nCells();

        // Gather/scatter list
        Pstream::gatherList(nCellsPerProc);
        Pstream::scatterList(nCellsPerProc);

        // Now that all processors have all the data, find the one with least
        // amount of cells (procMinI) and one with most cells (procMaxI)
        const label procMinI = findMin(nCellsPerProc);
        const label procMaxI = findMax(nCellsPerProc);

        // Convert number of cells to scalar
        const scalar minCells = nCellsPerProc[procMinI];
        const scalar maxCells = nCellsPerProc[procMaxI];

        // Calculate the imbalance
        const scalar imbalance = minCells/maxCells;

        // Print out the information
        Info<< "Minimum number of cells per processor: " << minCells << nl
            << "Maximum number of cells per processor: " << maxCells << nl
            << "Imbalance: " << imbalance << endl;

        // Perform load balancing only if the imbalance exceeds the threshold
        if (imbalance < imbalanceThreshold_)
        {
            Info<< "Performing load balancing..." << endl;

            // Propagate basic controls from refinementDict to loadBalanceDict,
            // hard coding the number of subdomains to be the same as number of
            // processors used in this run
            dictionary loadBalanceDict;

            // Get decomposition method and add it to the load balancing dict
            const word decompMethod = refinementDict_.lookup("method");
            loadBalanceDict.add(word("method"), decompMethod);

            // Set numberOfSubdomains to number or processors we are running
            loadBalanceDict.add(word("numberOfSubdomains"), Pstream::nProcs());

            // Perform load balancing
            const bool balanced = loadBalance(loadBalanceDict);

            hasChanged = hasChanged || balanced;
        }
        else
        {
            Info<< "No need to perform load balancing..." << endl;
        }
    }

    // Execute dummy mesh motion for the background mesh
    const pointField oldPoints = allPoints();
    fvMesh::movePoints(oldPoints);

    return hasChanged;
}


// ************************************************************************* //
