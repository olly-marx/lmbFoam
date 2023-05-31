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

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "numericsFunctions.H"
#include "interfaceFvData.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::numericsFunctions::limitVelocityField
(
    volVectorField& UField,
    const dictionary& dict
)
{
    // Limit the velocity only if we find the limit entry in the dictionary
    if (dict.found("limitMag" + UField.name()))
    {
        const scalar limitMagU =
            readScalar(dict.lookup("limitMag" + UField.name()));

        const volScalarField magUField(mag(UField));
        const scalar maxMagU = max(magUField).value();

        Info<< "mag(" << UField.name() << "): max: "
            << returnReduce(maxMagU, maxOp<scalar>())
            << ", average: "
            << returnReduce
               (
                   magUField.weightedAverage(UField.mesh().V()).value(),
                   sumOp<scalar>()
               )/Pstream::nProcs();

        if (maxMagU > limitMagU)
        {
            // Limit the field
            UField.internalField() *=
                neg(magUField.internalField() - limitMagU)
              + pos(magUField.internalField() - limitMagU)*
                limitMagU/(magUField.internalField() + SMALL);
            UField.correctBoundaryConditions();
            Info << " ...limiting" << endl;
        }
        else
        {
            Info << endl;
        }
    }
}


void Foam::numericsFunctions::extrapolateInterfaceVelocity
(
    volVectorField& UField,
    const dictionary& dict
)
{
    // Check whether to extrapolate the velocity field. Off by default
    const Switch extrapolateVelocity =
        dict.lookupOrDefault("extrapolateVelocity", false);

    if (extrapolateVelocity)
    {
        Info<< "Extrapolating field " << UField.name()
            << " near the interface." << endl;

        // Get necessary mesh data
        const fvMesh& mesh = UField.mesh();
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // Get interface data object
        const interfaceFvData& intFvData =
            mesh.lookupObject<interfaceFvData>("interfaceFvData");

        // Get a list of interface faces
        const interfaceFvData::DynamicLabelList& interfaceFaces =
            intFvData.interfaceFaces();

        // Get wet owners (1 = wet, 0 = dry)
        const surfaceScalarField& wetOwners = intFvData.wetOwners();

        // Get internal velocity field
        vectorField& UFieldIn = UField.internalField();

        // Prepare for extrapolation: count number of interface faces for each
        // cell and set velocity field to zero where needed
        labelField nInterfaceFacesPerCell(mesh.nCells(), 0);
        forAll (interfaceFaces, ifI)
        {
            // Get face index
            const label& faceI = interfaceFaces[ifI];

            // Get wet owner of the face
            const scalar& wetOwn = wetOwners[faceI];

            if (equal(wetOwn, 1.0))
            {
                // Get neighbour index
                const label& N = nei[faceI];

                // Owner of the face is wet, reset neighbour velocity
                UFieldIn[N] = vector::zero;
                
                // Increment number of interface faces for the dry cell
                ++nInterfaceFacesPerCell[N];
            }
            else
            {
                // Get owner index
                const label& P = own[faceI];

                // Neighbour of the face is wet, reset owner velocity
                UFieldIn[P] = vector::zero;

                // Increment number of interface faces for the dry cell
                ++nInterfaceFacesPerCell[P];
            }
        }

        // Loop through interface faces and correct the velocity from dry
        // cell to be average of the velocities from surrounding wet cells
        forAll (interfaceFaces, ifI)
        {
            // Get face index
            const label& faceI = interfaceFaces[ifI];

            // Get wet owner of the face
            const scalar& wetOwn = wetOwners[faceI];

            // Get owner and neighbour labels
            const label& P = own[faceI];
            const label& N = nei[faceI];

            if (equal(wetOwn, 1.0))
            {
                // Owner of the face is wet, set neighbour from owner value
                UFieldIn[N] += UFieldIn[P]/nInterfaceFacesPerCell[N];
            }
            else
            {
                // Neighbour of the face is wet, set owner from neighbour value
                UFieldIn[P] += UFieldIn[N]/nInterfaceFacesPerCell[P];
            }
        }

        // Switch to control second order extrapolation correction
        const Switch secondOrderExtrapolation =
            dict.lookupOrDefault("secondOrderExtrapolation", true);

        if (secondOrderExtrapolation)
        {
            // Now that the velocities are first-order extrapolated, we
            // calculate the second order update based on the gradient of
            // velocity field
            const volTensorField gradU = fvc::grad(UField);
            const tensorField& gradUIn = gradU.internalField();

            // Get cell centres
            const vectorField& CIn = mesh.C().internalField();

            // Now that the gradient is calculated, reset velocity field again
            // to allow average calculation
            forAll (interfaceFaces, ifI)
            {
                // Get face index
                const label& faceI = interfaceFaces[ifI];

                // Get wet owner of the face
                const scalar& wetOwn = wetOwners[faceI];

                if (equal(wetOwn, 1.0))
                {
                    // Owner of the face is wet, reset neighbour velocity
                    UFieldIn[nei[faceI]] = vector::zero;
                }
                else
                {
                    // Neighbour of the face is wet, reset owner velocity
                    UFieldIn[own[faceI]] = vector::zero;
                }
            }

            // Loop through interface faces
            forAll (interfaceFaces, ifI)
            {
                // Get face index
                const label& faceI = interfaceFaces[ifI];

                // Get wet owner of the face
                const scalar& wetOwn = wetOwners[faceI];

                // Get owner and neighbour
                const label& P = own[faceI];
                const label& N = nei[faceI];

                if (equal(wetOwn, 1.0))
                {
                    // Wet owner, extrapolate from owner to neighbour
                    UFieldIn[N] +=
                        (UFieldIn[P] + ((CIn[N] - CIn[P]) & gradUIn[P]))/
                        nInterfaceFacesPerCell[N];
                }
                else
                {
                    // Wet neighbour, extrapolate from neighbour to owner
                    UFieldIn[P] +=
                        (UFieldIn[N] + ((CIn[P] - CIn[N]) & gradUIn[N]))/
                        nInterfaceFacesPerCell[P];
                }
            }
        }

        // Correct boundary conditions
        UField.correctBoundaryConditions();

        // Note: This correction is not done across processor boundaries because
        // this function is currently in the testing phase and I assume that it
        // will be rarely used. If someone indeed starts to use this on a
        // regular basis, then I will parallelise it properly. VV, 17/Apr/2018.
    }
}


// ************************************************************************* //
