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

\*---------------------------------------------------------------------------*/

#include "interfaceFacesWave.H"
#include "polyMesh.H"
#include "wallPoint.H"
#include "MeshWave.H"
#include "globalMeshData.H"
#include "interfaceFvData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and number of iterations
Foam::interfaceFacesWave::interfaceFacesWave
(
    const polyMesh& mesh
)
:
    cellDistFuncs(mesh),
    distance_(mesh.nCells())
{
    // Get interfaceFvData object
    const interfaceFvData& intFvData =
        mesh.lookupObject<interfaceFvData>("interfaceFvData");

    // Get list of internal interface faces and patch interface faces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFaces();
    const interfaceFvData::DynamicLabelListList& bInterfaceFaces =
        intFvData.boundaryInterfaceFaces();

    // Get total number of faces
    label nFaces = interfaceFaces.size();
    forAll (bInterfaceFaces, patchI)
    {
        nFaces += bInterfaceFaces[patchI].size();
    }

    // Create lists which will hold the data. Note: using wallPoint to calculate
    // the distance
    List<wallPoint> faceDistances(nFaces);
    labelList changedFaces(nFaces);

    // Set changed faces and count them
    label nChangedFaces = 0;

    // Get interface points for interface faces
    const surfaceVectorField& interfacePoints = intFvData.cGamma();
    const vectorField& interfacePointsIn = interfacePoints.internalField();

    // Get mesh face centres (from fvMesh)
    const surfaceVectorField& faceCentres = intFvData.mesh().Cf();
    const vectorField& faceCentresIn = faceCentres.internalField();

    // Loop through internal faces first
    forAll (interfaceFaces, i)
    {
        // Get face index
        const label& meshFaceI = interfaceFaces[i];

        // Mark the face
        changedFaces[nChangedFaces] = meshFaceI;

        // Initialise the position of the neareast centre and squared distance
        // for this face
        faceDistances[nChangedFaces] =
            wallPoint
            (
                interfacePointsIn[meshFaceI],
                magSqr(interfacePointsIn[meshFaceI] - faceCentresIn[meshFaceI])
            );

        // Increment number of changed faces
        ++nChangedFaces;
    }

    // Get boundary fields
    const surfaceVectorField::GeometricBoundaryField& bInterfacePoints =
        interfacePoints.boundaryField();
    const surfaceVectorField::GeometricBoundaryField& bFaceCentres =
        faceCentres.boundaryField();

    // Loop through boundary faces
    forAll(bInterfaceFaces, patchI)
    {
        // Get interface faces on this patch
        const interfaceFvData::DynamicLabelList& pIntFaces =
            bInterfaceFaces[patchI];

        // Get polyPatch
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        // Get reference to patch fields
        const vectorField& pIntPoints = bInterfacePoints[patchI];
        const vectorField& pFaceCentres = bFaceCentres[patchI];

        // Loop through patch faces
        forAll (pIntFaces, i)
        {
            // Get patch local face index
            const label pfaceI = pIntFaces[i];

            // Mark the face with global index
            changedFaces[nChangedFaces] = patch.start() + pfaceI;

            // Initialise the position of the nearest centre and squared
            // distance for this face
            faceDistances[nChangedFaces] =
                wallPoint
                (
                    pIntPoints[pfaceI],
                    magSqr(pIntPoints[pfaceI] - pFaceCentres[pfaceI])
                );

            // Increment number of changed faces
            ++nChangedFaces;
        }

    }

    // Now the data is complete:
    // 1. changedFaces holds all global face indices to start from
    // 2. faceDistances holds all information necessary for calulating the least
    //    distance

    // Create the MeshWave object which will actually calculate the distance
    MeshWave<wallPoint> waveInfo
    (
        mesh,
        changedFaces,
        faceDistances,
        mesh.globalData().nTotalCells() // max number of iterations
    );

    // Get the distances from the mesh wave
    const List<wallPoint>& cellInfo = waveInfo.allCellInfo();

    // Copy cell distances into the placeholder field
    distance_.setSize(cellInfo.size());

    forAll (cellInfo, cellI)
    {
        // Get the distance
        const scalar dist = cellInfo[cellI].distSqr();

        if (cellInfo[cellI].valid())
        {
            // Cell has been covered, set the distance
            distance_[cellI] = Foam::sqrt(dist);
        }
        else
        {
            // Cell hasn't been visited, set the distance to -1
            distance_[cellI] = -1.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceFacesWave::~interfaceFacesWave()
{}


// ************************************************************************* //
