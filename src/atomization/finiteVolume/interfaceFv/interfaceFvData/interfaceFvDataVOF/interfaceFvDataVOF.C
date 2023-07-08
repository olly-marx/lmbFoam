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

\*---------------------------------------------------------------------------*/

#include "interfaceFvDataVOF.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceFvDataVOF::interfaceFvDataVOF
(
    const fvMesh& mesh,
    const volScalarField& alpha
)
:
    interfaceFvData(mesh),
    alpha_(alpha)
{
    // Update on construction
    this->update();

    // Update old data as well
    this->updateOld();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceFvDataVOF::update()
{
    // Bugfix: coupled patches need to be evaluated before patchNeighbourField
    // call (for interface faces on coupled patches). Rather use const_cast than
    // make alpha_ member a non-const reference. VV, 4/Mar/2015.
    alpha_.boundaryField().updateCoupledPatchFields();

    // First update old data
    this->updateOld();

    // Get references to const mesh data
    const fvMesh& mesh = this->mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get internal alpha field
    const scalarField& alphaIn = alpha_.internalField();

    // First update wetCells field
    this->wetCells() = pos(alphaIn - 0.5);

    // Get wet owners
    surfaceScalarField& wetOwn = this->wetOwners();
    scalarField& wetOwnIn = wetOwn.internalField();
    // Get wet neighbours
    surfaceScalarField& wetNei = this->wetNeighbours();
    scalarField& wetNeiIn = wetNei.internalField();

    // Get reference to interface faces
    DynamicLabelList& intFaces = this->interfaceFaces();
    DynamicLabelListList& boundaryIntFaces = this->boundaryInterfaceFaces();

    // Reset interface faces. Do not clear the storage to prevent excessive
    // resizing.
    intFaces.clear();
    forAll(boundaryIntFaces, patchI)
    {
        boundaryIntFaces[patchI].clear();
    }

    // Get reference to cell centres
    const volVectorField& C = this->mesh().C();
    const vectorField& CIn = C.internalField();

    // Calculate lambdas and cg at the faces before calculation of kappa (some
    // curvature models, e.g. interfaceFacesLevelSet rely on cg to be
    // up-to-date)
    surfaceVectorField& cg = this->cGamma();
    vectorField& cgIn = cg.internalField();
    surfaceScalarField lambdas
    (
        IOobject
        (
            "lambdas",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("minusOne", dimless, -1)
    );
    scalarField& lambdasIn = lambdas.internalField();

    // Populate internal wetOwners and interface faces
    forAll(neighbour, faceI)
    {
        // Get owner/neighbour of the face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get owner/neighbour alpha values
        const scalar& alphaOwn = alphaIn[own];
        const scalar& alphaNei = alphaIn[nei];

        // Populate wetOwners
        wetOwnIn[faceI] = pos(alphaOwn - 0.5);
        wetNeiIn[faceI] = pos(alphaNei - 0.5);

        // Level set changes sign across the interface face.
        if ((alphaOwn - 0.5)*(alphaNei - 0.5) < 0)
        {
            // Insert face index into dynamic list
            intFaces.append(faceI);

            // Calculate lambda for this face
            scalar& lambdaf = lambdasIn[faceI];
            lambdaf = (alphaOwn - 0.5)/(alphaOwn - alphaNei);

            // Calculate interface location for this face
            cgIn[faceI] = lambdaf*CIn[nei] + (1.0 - lambdaf)*CIn[own];
        }
    }

    // Populate boundary wetOwners and interface faces
    forAll(wetOwn.boundaryField(), patchI)
    {
        // Get reference to patch data
        const fvPatchScalarField& palpha = alpha_.boundaryField()[patchI];
        const fvPatch& patch = palpha.patch();

        fvsPatchScalarField& pwetOwn = wetOwn.boundaryField()[patchI];
        fvsPatchVectorField& pcg = cg.boundaryField()[patchI];
        fvsPatchScalarField& plambdas = lambdas.boundaryField()[patchI];
        DynamicLabelList& pIntFaces = boundaryIntFaces[patchI];

        // Get reference to face cells
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        // Populate wetOwners at the patch irrespective of the patch type
        forAll(mesh.boundary()[patchI], pfaceI)
        {
            pwetOwn[pfaceI] = pos(alphaIn[pFaceCells[pfaceI]] - 0.5);
        }

        // Populate interface faces at coupled patches
        if (palpha.coupled())
        {
            // Fetch alpha neighbour field and delta vectors
            const scalarField palphaNei = palpha.patchNeighbourField();
            const vectorField pd = patch.delta();

            forAll(palpha, pfaceI)
            {
                // Get alpha values
                const label& pcellI = pFaceCells[pfaceI];
                const scalar& alphaOwn = alphaIn[pcellI];
                const scalar& alphaNei = palphaNei[pfaceI];

                if ((alphaOwn - 0.5)*(alphaNei - 0.5) < 0)
                {
                    // Append coupled patch interface face
                    pIntFaces.append(pfaceI);

                    // Calculate lambda for this face
                    scalar& plambdaf = plambdas[pfaceI];
                    plambdaf = (alphaOwn - 0.5)/(alphaOwn - alphaNei);

                    // Calculate interface location for this face
                    pcg[pfaceI] = CIn[pcellI] + plambdaf*pd[pfaceI];
                }
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015
    }

    // Print number of internal faces
    Info<< "Number of internal interface faces: "
        << returnReduce(intFaces.size(), sumOp<label>()) << endl;

    // Grab blending fields at interface faces
    surfaceScalarField& bOverbar = this->betaOverbar();
    scalarField& bOverbarIn = bOverbar.internalField();
    surfaceScalarField& hJump = this->hydrostaticJump();
    scalarField& hJumpIn = hJump.internalField();

    // bOverbar and hJump reset is ommited - the values for faces that will not
    // be overwritten are not used

    // Get betaPlus and betaMinus values
    const scalar bPlus = this->betaPlus();
    const scalar bMinus = this->betaMinus();

    // Calculate the curvature and get internal field
    const volScalarField kappaInterface = curvature().kappa(alpha_);
    // Update coupled boundaries.  Does this belong somewhere else?
    kappaInterface.boundaryField().updateCoupledPatchFields();
    
    const scalarField& kappaIn = kappaInterface.internalField();

    // Calculate betaOverbar and hydrostaticJump at internal interface faces
    forAll (intFaces, ifI)
    {
        // Get face index
        const label& faceI = intFaces[ifI];

        // Get owner/neighbour labels for cells sharing this face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get lambda for this face
        const scalar& lambdaf = lambdasIn[faceI];

        // Get wet owner blending field (1 = wet, 0 = dry)
        const scalar& wetOwner = wetOwnIn[faceI];

        // Calculate betaOverbar for this face
        bOverbarIn[faceI] =
            wetOwner*(bMinus*lambdaf + bPlus*(1 - lambdaf))
          + (1 - wetOwner)*(bPlus*lambdaf + bMinus*(1 - lambdaf));

        // Calculate curvature at the interface location
        const scalar kappaf =
            lambdaf*kappaIn[nei] + (1 - lambdaf)*kappaIn[own];

        // Calculate jump term on this face
        hJumpIn[faceI] =
            // Jump due to density differences in gravitational field
            (1/bPlus - 1/bMinus)*(this->g() & cgIn[faceI])
            // Jump due to surface tension effects
          - this->sigma().value()*kappaf;
    }

    // Calculate betaOverbar and hydrostaticJump at boundary interface faces
    forAll (boundaryIntFaces, patchI)
    {
        // Get references to patch data
        const fvPatchScalarField& palpha = alpha_.boundaryField()[patchI];
        const fvPatchScalarField& pkappa =
            kappaInterface.boundaryField()[patchI];
        const fvPatch& patch = palpha.patch();
        const vectorField pd = patch.delta();
        const fvsPatchScalarField& pwetOwn = wetOwn.boundaryField()[patchI];
        const DynamicLabelList& pIntFaces = boundaryIntFaces[patchI];
        const fvsPatchScalarField& plambdas = lambdas.boundaryField()[patchI];
        const fvsPatchVectorField& pcg = cg.boundaryField()[patchI];

        fvsPatchScalarField& pbOverbar = bOverbar.boundaryField()[patchI];
        fvsPatchScalarField& phJump = hJump.boundaryField()[patchI];

        // Get reference to face cells
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        // Calculation dependant on patch type
        if (palpha.coupled())
        {
            // Fetch alpha and kappa neighbour field
            const scalarField palphaNei = palpha.patchNeighbourField();
            const scalarField pkappaNei = pkappa.patchNeighbourField();

            forAll(pIntFaces, pifI)
            {
                // Get face and cell index
                const label& pfaceI = pIntFaces[pifI];
                const label& pcellI = pFaceCells[pfaceI];

                // Get lambda for this face
                const scalar plambdaf = plambdas[pfaceI];

                // Get owner blending field for this patch face
                const scalar& pwetOwner = pwetOwn[pfaceI];

                // Calculate betaOverbar for this patch face
                pbOverbar[pfaceI] =
                    pwetOwner*(bMinus*plambdaf + bPlus*(1 - plambdaf))
                  + (1 - pwetOwner)*(bPlus*plambdaf + bMinus*(1 - plambdaf));

                // Calculate curvature at the interface location
                const scalar pkappaf =
                    plambdaf*pkappaNei[pfaceI] + (1 - plambdaf)*kappaIn[pcellI];

                // Calculate jump term on this patch face
                phJump[pfaceI] =
                    // Jump due to density differences in gravitational field
                    (1/bPlus - 1/bMinus)*(this->g() & pcg[pfaceI])
                    // Jump due to surface tension effects
                  - this->sigma().value()*pkappaf;
            }
        }
        else
        {
            // MISSING CONTACT ANGLE EFFECT!  HJ, 30/Apr/2021
            // WIP!!!  This is still not working
            phJump = - this->sigma().value()*pkappa.snGrad();
            // forAll (pFaceCells, pfaceI)
            // {
            //     const label& pcellI = pFaceCells[pfaceI];

            //     // Get lambda for this face
            //     const scalar plambdaf = plambdas[pfaceI];

            //     // Calculate curvature at the interface location
            //     const scalar pkappaf =
            //         plambdaf*pkappa[pfaceI] + (1 - plambdaf)*kappaIn[pcellI];

            //     phJump[pfaceI] =
            //         // Jump due to surface tension effects
            //         - this->sigma().value()*pkappaf;
            // }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015

    }
}


// ************************************************************************* //
