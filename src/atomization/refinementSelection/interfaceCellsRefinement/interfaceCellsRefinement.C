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

#include "interfaceCellsRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "interfaceFvData.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interfaceCellsRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    interfaceCellsRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCellsRefinement::interfaceCellsRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    nLayers_(readLabel(coeffDict().lookup("nLayers"))),
    pointBasedExtension_
    (
        coeffDict().lookupOrDefault<Switch>("pointBasedExtension", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::interfaceCellsRefinement::~interfaceCellsRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::interfaceCellsRefinement::refinementCellCandidates() const
{
    // Get interface data
    const interfaceFvData& intFvData =
        mesh().lookupObject<interfaceFvData>("interfaceFvData");

    // Mark-up for cells sharing interface faces
    boolList refCells(mesh().nCells(), false);

    // Get mesh data
    const unallocLabelList& nei = mesh().neighbour();
    const unallocLabelList& own = mesh().owner();

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFaces();

    // Collect cells near internal interface faces
    forAll (interfaceFaces, i)
    {
        // Get face index
        const label& faceI = interfaceFaces[i];

        // Mark owner and neighbour cells
        refCells[own[faceI]] = true;
        refCells[nei[faceI]] = true;
    }

    // Get reference to the list of boundary interface faces
    const interfaceFvData::DynamicLabelListList& bInterfaceFaces =
        intFvData.boundaryInterfaceFaces();

    const fvBoundaryMesh& bMesh = mesh().boundary();

    // Loop through all patches
    forAll (bInterfaceFaces, patchI)
    {
        const fvPatch& patch = bMesh[patchI];

        // Only coupled patches needed
        if (patch.coupled())
        {
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                bInterfaceFaces[patchI];

            const unallocLabelList& pFaceCells = patch.faceCells();

            // Loop through all interface faces on this patch
            forAll (pinterfaceFaces, pifI)
            {
                const label& pfaceI = pinterfaceFaces[pifI];

                refCells[pFaceCells[pfaceI]] = true;
            }
        }
    }

    // Extend marked cells nLayers
    for (label i = 0; i < nLayers_; ++i)
    {
        if (pointBasedExtension_)
        {
            // Note: does syncying across processor points
            meshTools::extendMarkedCellsAcrossPoints(mesh(), refCells);
        }
        else
        {
            // Note: does syncying across processor faces
            meshTools::extendMarkedCellsAcrossFaces(mesh(), refCells);
        }
    }

    // Collect all cells to refine
    dynamicLabelList refinementCandidates(mesh().nCells()/100);
    forAll (refCells, cellI)
    {
        if (refCells[cellI])
        {
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
Foam::interfaceCellsRefinement::unrefinementPointCandidates() const
{
    // Get refinement cell candidates
    const labelList refinementCandidates = this->refinementCellCandidates();

    // Create point mask that will eventually be used to mark unrefinement
    // candidates
    boolList unrefinementMask(mesh().nPoints(), true);

    // Remove all points of all cells that are selected for refinement
    const labelListList& meshCellPoints = mesh().cellPoints();

    // Loop through all refinement cells
    forAll (refinementCandidates, i)
    {
        // Get cell index
        const label& cellI = refinementCandidates[i];

        // Loop through all points of this cell and mark them
        const labelList& cPoints = meshCellPoints[cellI];
        forAll (cPoints, cpI)
        {
            unrefinementMask[cPoints[cpI]] = false;
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

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //
