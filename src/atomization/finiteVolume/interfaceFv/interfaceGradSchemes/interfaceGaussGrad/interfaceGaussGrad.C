/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "interfaceGaussGrad.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceGaussGrad, 0);
gradScheme<scalar>::addIstreamConstructorToTable<interfaceGaussGrad>
    addInterfaceGaussGradScalarIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> interfaceGaussGrad::calcGrad
(
    const volScalarField& pd,
    const word& name
) const
{
    // Get reference to mesh
    const fvMesh& mesh = this->mesh();

    // Get reference to intFvData object
    const interfaceFvData& intFvData =
        mesh.lookupObject<interfaceFvData>("interfaceFvData");

    // Construct volVector gradient field
    tmp<volVectorField> tgGrad
    (
        new volVectorField
        (
            IOobject
            (
                "grad(" + pd.name() + ')',
                pd.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "0",
                pd.dimensions()/dimLength/intFvData.rhoPlus().dimensions(),
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        )
    );
    volVectorField& gGrad = tgGrad();

    // Get owner/neighbour addressing, surface area vectors and cell volumes
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const vectorField& SfIn = Sf.internalField();
    const scalarField& VIn = mesh.V().field();

    // Get weight field
    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(pd);
    const scalarField& wIn = tweights().internalField();

    // Get internal fields of pd and the gradient
    const scalarField& pdIn = pd.internalField();
    vectorField& igGrad = gGrad.internalField();

    // Ordinary Gauss grad for internal faces
    forAll(neighbour, faceI)
    {
        // Get owner/neighbour labels
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get weight at this face
        const scalar& wf = wIn[faceI];

        // Calculate Gauss grad contribution for this face
        const vector Sfpdf = SfIn[faceI]*(wf*pdIn[own] + (1 - wf)*pdIn[nei]);

        // Add owner/neighbour contributions
        igGrad[own] += Sfpdf;
        igGrad[nei] -= Sfpdf;
    }

    // Bugfix: Interpolate pd for correct coupled boundary treatment using the
    // same interpolation scheme as the one for gradient weights.
    // IG, 15/Jan/2018
    const surfaceScalarField interpolatedPd =
        this->tinterpScheme_().interpolate(pd);

    // Ordinary Gauss grad for boundary faces
    forAll(mesh.boundary(), patchI)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchI];
        // Bugfix: use interpolated pd for correct coupled boundary treatment,
        // IG, 15/Jan/2018
        const scalarField& ppd = interpolatedPd.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            igGrad[pFaceCells[faceI]] += pSf[faceI]*ppd[faceI];
        }
    }

    // STAGE 1: eliminate contributions across the interface from Gauss grad

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFacesOld();

    // Internal interface faces
    forAll(interfaceFaces, ifI)
    {
        // Get face index from list of interface faces
        const label& faceI = interfaceFaces[ifI];

        // Get owner/neighbour labels for cells sharing this face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get weight at this face
        const scalar& wf = wIn[faceI];

        // Calculate Gauss grad contribution for this face
        const vector Sfpdf = SfIn[faceI]*(wf*pdIn[own] + (1 - wf)*pdIn[nei]);

        // Subtract these contributions from owner and neighbour
        igGrad[own] -= Sfpdf;
        igGrad[nei] += Sfpdf;
    }

    // Also for coupled boundaries
    const interfaceFvData::DynamicLabelListList& boundaryInterfaceFaces =
        intFvData.boundaryInterfaceFacesOld();

    // Coupled patch interface faces
    forAll(boundaryInterfaceFaces, patchI)
    {
        if (pd.boundaryField()[patchI].coupled())
        {
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                boundaryInterfaceFaces[patchI];

            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();
            const vectorField& pSf = Sf.boundaryField()[patchI];
            // Bugfix: use interpolated pd for correct coupled boundary
            // treatment consistent with bugfix above. VV, 27/Jan/2019
            const scalarField& ppd = interpolatedPd.boundaryField()[patchI];

            forAll(pinterfaceFaces, pifI)
            {
                // Get face index from list of interface faces
                const label& pfaceI = pinterfaceFaces[pifI];

                // Subract the contribution to face cell
                igGrad[pFaceCells[pfaceI]] -= pSf[pfaceI]*ppd[pfaceI];
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015.
    }

    // STAGE 2: density (1/rho) jump condition

    // Get necessary constants from interfaceFvData object
    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Get wet cells field
    const scalarField& wetCells = intFvData.wetCellsOld();

    // Scale the gradient in cell with appropriate rho - internal field
    // contributions
    igGrad *= wetCells*betaPlus + (1 - wetCells)*betaMinus;

    // Boundary contributions not needed because it is the cell centred field

    // STAGE 3: pressure gradient jump condition

    // Get references to additional interface finite volume data
    const surfaceScalarField& betaOverbar = intFvData.betaOverbarOld();
    const scalarField& betaOverbarIn = betaOverbar.internalField();
    const surfaceScalarField& hJump = intFvData.hydrostaticJumpOld();
    const scalarField& hJumpIn = hJump.internalField();

    // Get reference to wet owners (wet owner = 1, dry owner = 0)
    const surfaceScalarField& wetOwners = intFvData.wetOwnersOld();
    const scalarField& wetOwnersIn = wetOwners.internalField();

    // Add jump conditions due to internal interface faces
    forAll(interfaceFaces, ifI)
    {
        // Get face index from list of interface faces
        const label& faceI = interfaceFaces[ifI];

        // Get owner/neighbour labels for cells sharing this face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get ordinary weight at this face
        const scalar& wf = wIn[faceI];

        // Get owner/neighbour pd, cell volumes and this surface area vector
        const scalar& pdOwn = pdIn[own];
        const scalar& pdNei = pdIn[nei];
        const vector& SfInf = SfIn[faceI];

        // Helper variable: betaPlus*betaMinus/betaOverbar
        const scalar betaBlend = betaPlus*betaMinus/betaOverbarIn[faceI];

        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Helper: delta pd with jump conditions
        const scalar deltapd = pdNei - pdOwn - (2*wetOwner - 1)*hJumpIn[faceI];

        // Add owner contribution
        igGrad[own] += SfInf*
        (
            // First order extrapolation
            (wetOwner*betaPlus + (1 - wetOwner)*betaMinus)*pdOwn
            // Second order correction
          + betaBlend*(1 - wf)*deltapd
        );

        // Add neighbour contribution
        igGrad[nei] -= SfInf*
        (
            // First order extrapolation
            (wetOwner*betaMinus + (1 - wetOwner)*betaPlus)*pdNei
            // Second order correction
          - betaBlend*wf*deltapd
        );
    }

    // Add jump conditions due to boundary interface faces
    forAll(boundaryInterfaceFaces, patchI)
    {
        const fvPatch& patch = wetOwners.boundaryField()[patchI].patch();

        if (patch.coupled())
        {
            // Get references to interface related patch data
            const fvsPatchScalarField& pwetOwners =
                wetOwners.boundaryField()[patchI];
            const fvsPatchScalarField& pbetaOverbar =
                betaOverbar.boundaryField()[patchI];
            const fvsPatchScalarField& phJump = hJump.boundaryField()[patchI];
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                boundaryInterfaceFaces[patchI];

            // Get references to additional fields
            const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];
            const scalarField ppdNeiField =
                pd.boundaryField()[patchI].patchNeighbourField();
            const vectorField& pSf = Sf.boundaryField()[patchI];

            // Get reference to face cells
            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(pinterfaceFaces, pifI)
            {
                // Get face and cell index
                const label& pfaceI = pinterfaceFaces[pifI];
                const label& pcellI = pFaceCells[pfaceI];

                // Get weight at this face
                const scalar& pwf = pw[pfaceI];

                // Get pd values, face cell volume and surface area vector
                const scalar& ppdOwn = pdIn[pcellI];
                const scalar& ppdNei = ppdNeiField[pfaceI];
                const vector& pSff = pSf[pfaceI];

                // Helper variable: betaPlus*betaMinus/betaOverbar
                const scalar pbetaBlend =
                    betaPlus*betaMinus/pbetaOverbar[pfaceI];

                // Get wet owner blending
                const scalar& pwetOwner = pwetOwners[pfaceI];

                // Helper: delta pd with jump conditions
                const scalar pdeltapd =
                    ppdNei - ppdOwn - (2*pwetOwner - 1)*phJump[pfaceI];

                // Add owner contribution
                igGrad[pcellI] += pSff*
                (
                    // First order extrapolation
                    (pwetOwner*betaPlus + (1 - pwetOwner)*betaMinus)*ppdOwn
                    // Second order correction
                  + pbetaBlend*(1 - pwf)*pdeltapd
                );
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015
    }

    // Divide with cell volume
    igGrad /= VIn;

    // Correct boundary conditions - first reset using zeroGradient and then
    // correct using snGrad of pd boundary field
    gGrad.correctBoundaryConditions();
    interfaceGradScheme::interfaceCorrectBoundaryConditions
    (
        pd,
        gGrad,
        intFvData
    );

    return tgGrad;
}





tmp<BlockLduSystem<vector, vector> >
interfaceGaussGrad::fvmGrad( const volScalarField& pd) const
{
    // Get reference to mesh
    const fvMesh& mesh = pd.mesh();

    // Get reference to intFvData object
    const interfaceFvData& intFvData =
        mesh.lookupObject<interfaceFvData>("interfaceFvData");

    // Get owner/neighbour addressing, surface area vectors and cell volumes
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const vectorField& SfIn = Sf.internalField();

    // Get weight field
    tmp<surfaceScalarField> tweights = this->tinterpScheme_().weights(pd);
    const scalarField& wIn = tweights().internalField();

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    // Ordinary Gauss grad for internal faces
    l = -wIn*SfIn;
    u = l + SfIn;

//    forAll(u, I)
//    {
//        d[owner[I]] += SfIn[I]*wIn[I];
//        d[neighbour[I]] -= SfIn[I]*(1 - wIn[I]);
//    }


    bs.negSumDiag();

    // Ordinary boundary contributions
    forAll(pd.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = pd.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& Sf = patch.Sf();
        const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));
        // Diag contribution
        forAll(pf, faceI)
        {
            d[fc[faceI]] += internalCoeffs[faceI]*Sf[faceI];
        }

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const vectorField pcl = -pw*Sf;
            const vectorField pcu = pcl + Sf;

            // Coupling  contributions
            pcoupleLower -= pcl;
            pcoupleUpper -= pcu;
        }
        else
        {
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Boundary contribution
            forAll(pf, faceI)
            {
                source[fc[faceI]] -= boundaryCoeffs[faceI]*Sf[faceI];
            }
        }
    }

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFaces();

    // Multiply coefficents with 1/rho

    // Get necessary constants from interfaceFvData object
    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Get reference to wet owners and neighbours (wet owner = 1, dry owner = 0)
    const surfaceScalarField& wetOwners = intFvData.wetOwners();
    const scalarField& wetOwnersIn = wetOwners.internalField();
    const surfaceScalarField& wetNeighbours = intFvData.wetNeighbours();
    const scalarField& wetNeighboursIn = wetNeighbours.internalField();

    // Get wet cells field
    const scalarField& wetCells = intFvData.wetCells();

    // Helper variables
    const scalarField betaOwn
    (
        wetOwnersIn*betaPlus + (1 - wetOwnersIn)*betaMinus
    );
    const scalarField betaNei
    (
        wetNeighboursIn*betaPlus + (1 - wetNeighboursIn)*betaMinus
    );
    const scalarField betaCell
    (
        wetCells*betaPlus + (1 - wetCells)*betaMinus
    );

    // STAGE 2, 1/rho jump
    l *= betaNei;
    u *= betaOwn;
    d *= betaCell;
    source *= betaCell;

    // STAGE 3: pressure gradient jump condition

    // Get references to additional interface finite volume data
    const surfaceScalarField& betaOverbar = intFvData.betaOverbar();
    const scalarField& betaOverbarIn = betaOverbar.internalField();
    const surfaceScalarField& hJump = intFvData.hydrostaticJump();
    const scalarField& hJumpIn = hJump.internalField();

    // Add jump conditions due to internal interface faces
    forAll(interfaceFaces, ifI)
    {
        // Get face index from list of interface faces
        const label& faceI = interfaceFaces[ifI];

        // Get owner/neighbour labels for cells sharing this face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Get ordinary weight at this face
        const scalar& wf = wIn[faceI];

        // Get this surface area vector
        const vector& SfInf = SfIn[faceI];

        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Helper variables
        const scalar betaOverbarf = betaOverbarIn[faceI];
        const scalar betaBlend = betaPlus*betaMinus/betaOverbarf;
        const scalar betaMinusOverbar = betaMinus/betaOverbarf;
        const scalar betaPlusOverbar = betaPlus/betaOverbarf;
        const scalar betaOwnf = betaOwn[faceI];
        const scalar betaNeif = betaNei[faceI];

        l[faceI] *= (1 - wetOwner)*betaMinusOverbar + wetOwner*betaPlusOverbar;
        u[faceI] *= wetOwner*betaMinusOverbar + (1 - wetOwner)*betaPlusOverbar;

        // Remove diagonal contribution across the interface faces
        d[own] -= SfInf*wf*betaOwnf;
        d[nei] -= - SfInf*(1 - wf)*betaNeif;

        // Add jump conditions to diagonal coeffs
        d[own] += SfInf*
        (
            // First order contribution
            betaOwnf
            // Second order correction
          - betaBlend*(1 - wf)
        );

        d[nei] += - SfInf*
        (
            // First order contribution
            betaNeif
            // Second order correction
          - betaBlend*wf
        );

        source[own] +=  SfInf*betaBlend*(1 - wf)*
                        (2*wetOwner - 1)*hJumpIn[faceI];
        source[nei] +=  SfInf*betaBlend*wf*
                        (2*wetOwner - 1)*hJumpIn[faceI];
    }


    // Add 1/rho conditions for coupled patches

    forAll(pd.boundaryField(), patchI)
    {
        const fvPatch& patch = pd.boundaryField()[patchI].patch();

        if (patch.coupled())
        {
            CoeffField<vector>::linearTypeField& pcU =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcL =
                bs.coupleLower()[patchI].asLinear();

            const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

            // Multiply upper and lower coeffs with beta in cell center. This
            // assumes that beta is equal in the cell across the coupled face,
            // which is not true only in the case of interfaceFace, which will
            // be taken into account later.
            forAll(patch, faceI)
            {
                pcU[faceI] *= betaCell[pFaceCells[faceI]];
                pcL[faceI] *= betaCell[pFaceCells[faceI]];
            }
        }
    }
    // Add jump conditions for interface coupled patches
    const interfaceFvData::DynamicLabelListList& boundaryInterfaceFaces =
    intFvData.boundaryInterfaceFaces();

    forAll(boundaryInterfaceFaces, patchI)
    {
        const fvPatch& patch = wetOwners.boundaryField()[patchI].patch();

        if (patch.coupled())
        {

            CoeffField<vector>::linearTypeField& pcU =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcL =
                bs.coupleLower()[patchI].asLinear();

            const fvPatchScalarField& pf = pd.boundaryField()[patchI];
            const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];

            const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Get references to interface related patch data
            const fvsPatchScalarField& pwetOwners =
                wetOwners.boundaryField()[patchI];
            const fvsPatchScalarField& pbetaOverbar =
                betaOverbar.boundaryField()[patchI];
            const fvsPatchScalarField& phJump = hJump.boundaryField()[patchI];
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                boundaryInterfaceFaces[patchI];

            // Get references to additional fields
            const vectorField& pSf = Sf.boundaryField()[patchI];

            // Get reference to face cells
            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(pinterfaceFaces, pifI)
            {
                // Get face and cell index
                const label& pfaceI = pinterfaceFaces[pifI];
                const label& pcellI = pFaceCells[pfaceI];

                // Get surface area vector
                const vector& pSff = pSf[pfaceI];

                // Get wet owner blending
                const scalar& pwetOwner = pwetOwners[pfaceI];

                // Helper variables
                const scalar pbetaBlend =
                betaPlus*betaMinus/pbetaOverbar[pfaceI];

                // Internal coeffs
                const scalar iCp = internalCoeffs[pfaceI];
                // Boundary coeffs
                const scalar bCp = boundaryCoeffs[pfaceI];

                // Get beta for the boundary cell
                const scalar betaCellp = betaCell[pcellI];

                // Remove ordinary 1/rho conditions to upper and lower coeffs and
                // add jump conditions
                pcL[pfaceI] *= pbetaBlend/betaCell[pcellI];
                pcU[pfaceI] *= pbetaBlend/betaCell[pcellI];

                // Remove ordinary diagonal contribution
                d[pcellI] -= iCp*pSff*betaCellp;

                // Add jump condition contribution to the diagonal
                d[pcellI] += pSff*
                (
                    // First order contribution
                    betaCellp
                    // Second order correction
                  - pbetaBlend*iCp
                );

                // Add jump condition contribution to the source
                source[pcellI] += pSff*pbetaBlend*bCp*
                                (2*pwetOwner - 1)*phJump[pfaceI];
            }
        }

    }


    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
