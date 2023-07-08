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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "interfaceProximityRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "labelIOField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interfaceProximityRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    interfaceProximityRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProximityRefinement::interfaceProximityRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    refinementDistanceFront_
    (
        readScalar(coeffDict().lookup("refinementDistanceFront"))
    ),
    refinementDistanceBack_
    (
        readScalar(coeffDict().lookup("refinementDistanceBack"))
    ),
    unrefinementDistanceFront_
    (
        readScalar(coeffDict().lookup("unrefinementDistanceFront"))
    ),
    unrefinementDistanceBack_
    (
        readScalar(coeffDict().lookup("unrefinementDistanceBack"))
    ),
    relativeDistance_
    (
        (coeffDict().lookupOrDefault<Switch>("relativeDistance", true))
    ),
    interfaceWaveInfo_()
{
    // Sanity checks
    if (refinementDistanceFront_ < SMALL)
    {
        FatalIOErrorIn
        (
            "interfaceProximityRefinement::interfaceProximityRefinement"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict"
            "\n)",
            coeffDict()
        )   << "Specified zero or negative refinementDistanceFront." << nl
            << "The entry must be positive."
            << abort(FatalIOError);
    }

    if (refinementDistanceBack_ < SMALL)
    {
        FatalIOErrorIn
        (
            "interfaceProximityRefinement::interfaceProximityRefinement"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict"
            "\n)",
            coeffDict()
        )   << "Specified zero or negative refinementDistanceBack." << nl
            << "The entry must be positive."
            << abort(FatalIOError);
    }

    if (unrefinementDistanceFront_ < SMALL)
    {
        FatalIOErrorIn
        (
            "interfaceProximityRefinement::interfaceProximityRefinement"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict"
            "\n)",
            coeffDict()
        )   << "Specified zero or negative unrefinementDistanceFront." << nl
            << "The entry must be positive."
            << abort(FatalIOError);
    }

    if (unrefinementDistanceBack_ < SMALL)
    {
        FatalIOErrorIn
        (
            "interfaceProximityRefinement::interfaceProximityRefinement"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict"
            "\n)",
            coeffDict()
        )   << "Specified zero or negative unrefinementDistanceBack." << nl
            << "The entry must be positive."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::interfaceProximityRefinement::~interfaceProximityRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::interfaceProximityRefinement::refinementCellCandidates() const
{
    // Calculate the relative distance based on current refinement level
    const labelIOField& cellLevel =
        mesh().lookupObject<labelIOField>("cellLevel");

    // Note negative sign for frontal refinement distance since we will use
    // signed distance to the interface (i.e. Level Set)
    scalarField rDistFront(mesh().nCells(), -refinementDistanceFront_);
    scalarField rDistBack(mesh().nCells(), refinementDistanceBack_);
    if (relativeDistance_)
    {
        forAll (rDistFront, cellI)
        {
            const scalar relativeFactor = 1.0/pow(2, cellLevel[cellI]);

            rDistFront[cellI] *= relativeFactor;
            rDistBack[cellI] *= relativeFactor;
        }
    }

    // Make sure that the data is cleared-up from previous time-step or
    // iteration and re-calculate. Basically we assume that the interface data
    // used by interfaceFacesWave has been updated between successive calls to
    // refinementCellCandidates()
    interfaceWaveInfo_.clear();
    interfaceWaveInfo_.set(new interfaceFacesWave(mesh()));

    // Get alpha field
    const volScalarField& alpha = mesh().lookupObject<volScalarField>("alpha1");
    const scalarField& alphaIn = alpha.internalField();

    // Get signed distance to the interface (positive when alpha > 0.5, negative
    // when alpha < 0.5)
    const scalarField interfaceDistance =
        sign(alphaIn - 0.5)*interfaceWaveInfo_->distance();

    // Create storage for collection of cells. Assume that 10% of all
    // cells will be marked to prevent excessive resizing.
    dynamicLabelList refinementCandidates(mesh().nCells()/10);

    // Loop through cells and collect refinement candidates
    forAll (interfaceDistance, cellI)
    {
        const scalar& psi = interfaceDistance[cellI];

        if
        (
            (psi < rDistBack[cellI])
         && (psi > rDistFront[cellI])
        )
        {
            // Cell is near the interface, append it for potential refinement
            refinementCandidates.append(cellI);
        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " cells as refinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::interfaceProximityRefinement::unrefinementPointCandidates() const
{
    // Calculate the relative distance based on current refinement level
    const labelIOField& cellLevel =
        mesh().lookupObject<labelIOField>("cellLevel");

    // Note negative sign for frontal refinement distance since we will use
    // signed distance to the interface (i.e. Level Set)
    scalarField rDistFront(mesh().nCells(), -unrefinementDistanceFront_);
    scalarField rDistBack(mesh().nCells(), unrefinementDistanceBack_);
    if (relativeDistance_)
    {
        forAll (rDistFront, cellI)
        {
            const scalar relativeFactor = 1.0/pow(2, cellLevel[cellI]);

            rDistFront[cellI] *= relativeFactor;
            rDistBack[cellI] *= relativeFactor;
        }
    }

    // Make sure that the distance to interface is calculated
    if (interfaceWaveInfo_.empty())
    {
        interfaceWaveInfo_.set(new interfaceFacesWave(mesh()));
    }

    // Get alpha field
    const volScalarField& alpha = mesh().lookupObject<volScalarField>("alpha1");
    const scalarField& alphaIn = alpha.internalField();

    // Get signed distance to the interface (positive when alpha > 0.5, negative
    // when alpha < 0.5)
    const scalarField interfaceDistance =
        sign(alphaIn - 0.5)*interfaceWaveInfo_->distance();

    // Create mark-up field for all points touched by unrefinement cell
    // candidates (cells that are far away from the interface)
    boolList unrefinementMask(mesh().nPoints(), false);

    // Get cell-point addressing
    const labelListList& meshCellPoints = mesh().cellPoints();

    // Loop through cells and mask unrefinement point candidates
    forAll (interfaceDistance, cellI)
    {
        const scalar& psi = interfaceDistance[cellI];

        if
        (
            (psi > rDistBack[cellI])
         || (psi < rDistFront[cellI])
        )
        {
            // This cell is far away from the interface, get its points
            const labelList& cPoints = meshCellPoints[cellI];

            // Mark points as possible unrefinement candidates
            forAll (cPoints, i)
            {
                unrefinementMask[cPoints[i]] = true;
            }
        }
    }

    // List for collecting unrefinement candidates. Assume all points are going
    // to get unrefinemed to prevent excessive resizing
    dynamicLabelList unrefinementCandidates(mesh().nPoints());

    // Collect marked points into the list
    forAll (unrefinementMask, pointI)
    {
        if (unrefinementMask[pointI])
        {
            unrefinementCandidates.append(pointI);
        }
    }

    // Clear out the data member for next iteration
    interfaceWaveInfo_.clear();

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //
