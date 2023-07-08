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

#include "interfaceLeastSquaresGrad.H"
#include "leastSquaresVectors.H"
#include "interfaceGaussGrad.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceLeastSquaresGrad, 0);
gradScheme<scalar>::addIstreamConstructorToTable<interfaceLeastSquaresGrad>
    addInterfaceLeastSquaresGradScalarIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> interfaceLeastSquaresGrad::calcGrad
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
    tmp<volVectorField> tlsGrad
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
    volVectorField& lsGrad = tlsGrad();

    // Get reference to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const vectorField& ownLsIn = ownLs.internalField();
    const surfaceVectorField& neiLs = lsv.nVectors();
    const vectorField& neiLsIn = neiLs.internalField();

    // Get owner/neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get internal fields of pd and the gradient
    const scalarField& pdIn = pd.internalField();
    vectorField& ilsGrad = lsGrad.internalField();

    // Ordinary least squares grad for internal faces
    forAll(neighbour, faceI)
    {
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        const scalar deltapd = pdIn[nei] - pdIn[own];

        ilsGrad[own] += ownLsIn[faceI]*deltapd;
        ilsGrad[nei] -= neiLsIn[faceI]*deltapd;
    }

    // Ordinary least squares grad for boundary faces
    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchVectorField& pOwnLs = ownLs.boundaryField()[patchI];

        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        if (pd.boundaryField()[patchI].coupled())
        {
            const scalarField neipd =
                pd.boundaryField()[patchI].patchNeighbourField();

            forAll(neipd, pfaceI)
            {
                lsGrad[pFaceCells[pfaceI]] +=
                    pOwnLs[pfaceI]*(neipd[pfaceI] - pdIn[pFaceCells[pfaceI]]);
            }
        }
        else
        {
            const fvPatchScalarField& ppd = pd.boundaryField()[patchI];

            forAll(ppd, pfaceI)
            {
                lsGrad[pFaceCells[pfaceI]] +=
                     pOwnLs[pfaceI]*(ppd[pfaceI] - pdIn[pFaceCells[pfaceI]]);
            }
        }
    }

    // STAGE 1: eliminate contributions across the interface from least squares
    // grad

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFacesOld();

    // Internal interface faces
    forAll(interfaceFaces, ifI)
    {
        // Get face index from list of interface faces
        const label& faceI = interfaceFaces[ifI];

        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        const scalar deltapd = pdIn[nei] - pdIn[own];

        // Subtract contributions from owner and neighbour
        ilsGrad[own] -= ownLsIn[faceI]*deltapd;
        ilsGrad[nei] += neiLsIn[faceI]*deltapd;
    }

    // Also for coupled boundaries
    const interfaceFvData::DynamicLabelListList& boundaryInterfaceFaces =
        intFvData.boundaryInterfaceFacesOld();

    // Coupled patch interface faces
    forAll(boundaryInterfaceFaces, patchI)
    {
        const fvPatch& patch = pd.boundaryField()[patchI].patch();

        if (patch.coupled())
        {
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                boundaryInterfaceFaces[patchI];
            const fvsPatchVectorField& pOwnLs = ownLs.boundaryField()[patchI];

            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();
            const scalarField neipd =
                pd.boundaryField()[patchI].patchNeighbourField();

            forAll(pinterfaceFaces, pifI)
            {
                // Get face index from list of interface faces
                const label& pfaceI = pinterfaceFaces[pifI];

                // Subtract the contribution to face cell
                lsGrad[pFaceCells[pfaceI]] -=
                    pOwnLs[pfaceI]*(neipd[pfaceI] - pdIn[pFaceCells[pfaceI]]);
            }
        }
    }

    // STAGE 2: density (1/rho) jump condition

    // Get necessary constants from interfaceFvData object
    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Get wet cells field
    const scalarField& wetCells = intFvData.wetCellsOld();

    // Scale the gradient in cell with appropriate rho - internal field
    // contribution
    ilsGrad *= wetCells*betaPlus + (1 - wetCells)*betaMinus;

    // Boundary contributions not needed because it is cell cetnred field

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
        // Get face index from list of interface face
        const label& faceI = interfaceFaces[ifI];

        // Get owner/neighbour labels for cells sharing this face
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        // Helper: delta pd with jump conditions
        const scalar blendedDeltapd = betaPlus*betaMinus/betaOverbarIn[faceI]*
            (pdIn[nei] - pdIn[own] - (2*wetOwnersIn[faceI] - 1)*hJumpIn[faceI]);

        // Add owner/neighbour contributions
        ilsGrad[own] += ownLsIn[faceI]*blendedDeltapd;
        ilsGrad[nei] -= neiLsIn[faceI]*blendedDeltapd;
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

            // Get reference to neighbour field
            const scalarField neipd =
                pd.boundaryField()[patchI].patchNeighbourField();

            // get reference to face cells
            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(pinterfaceFaces, pifI)
            {
                // Get face and cell index
                const label& pfaceI = pinterfaceFaces[pifI];
                const label& pcellI = pFaceCells[pfaceI];

                // Helper: delta pd with jump conditions
                const scalar pblendedDeltapd =
                    betaPlus*betaMinus/pbetaOverbar[pfaceI]*
                    (
                        neipd[pfaceI] - pdIn[pcellI]
                     - (2*pwetOwners[pfaceI] - 1)*phJump[pfaceI]
                    );

                // Add owner contribution
                ilsGrad[pcellI] += ownLsIn[pcellI]*pblendedDeltapd;
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015
    }

    // Correct boundary conditions
    lsGrad.correctBoundaryConditions();
    interfaceGradScheme::interfaceCorrectBoundaryConditions
    (
        pd,
        lsGrad,
        intFvData
    );

    return tlsGrad;
}


tmp<BlockLduSystem<vector, vector> >
interfaceLeastSquaresGrad::fvmGrad( const volScalarField& pd) const
{
    const fvMesh& mesh = pd.mesh();

    // Get reference to intFvData object
    const interfaceFvData& intFvData =
        mesh.lookupObject<interfaceFvData>("interfaceFvData");

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // This to be moved to FV mesh.  HJ, 15/Oct/2021
    volScalarField cellV
    (
        IOobject
        (
            "cellV",
            pd.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0),
        "zeroGradient"
    );
    cellV.internalField() = mesh.V();
    cellV.correctBoundaryConditions();
    const scalarField& cellVIn = cellV.internalField();

    const surfaceScalarField& w = mesh.weights();

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

    // Get references to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    //STAGE 1: Ordinary least squares gradient
    forAll(neighbour, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        u[faceI] = cellVIn[own]*ownLs[faceI];
        l[faceI] = cellVIn[nei]*neiLs[faceI];
    }

    // Boundary contributions
    forAll(pd.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = pd.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& pownLs = ownLs.boundaryField()[patchI];
        const fvsPatchScalarField& pw = w.boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        // Part of diagonal contribution irrespective of the patch type
        forAll(pf, faceI)
        {
            const label cellI = fc[faceI];
            d[cellI] -= cellVIn[cellI]*pownLs[faceI];
        }

        if (patch.coupled())
        {
            const vectorField& pneiLs = neiLs.boundaryField()[patchI];
            const scalarField cellVInNei =
                cellV.boundaryField()[patchI].patchNeighbourField();

            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            // Coupling  and diagonal contributions
            forAll(pf, faceI)
            {
                const vector upper = cellVIn[fc[faceI]]*pownLs[faceI];
                const vector lower = cellVInNei[faceI]*pneiLs[faceI];

                pcoupleUpper[faceI] -= upper;
                pcoupleLower[faceI] -= lower;
            }
        }
        else
        {
            const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Diagonal and source contributions depending on the patch type
            forAll(pf, faceI)
            {
                const label cellI = fc[faceI];
                d[cellI] += cellVIn[cellI]*pownLs[faceI]*internalCoeffs[faceI];
                source[cellI] -= cellVIn[cellI]*pownLs[faceI]*
                    boundaryCoeffs[faceI];
            }
        }
    }

    // STAGE 2: 1/rho jump condition
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

    // 1/rho jump
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

        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Helper variables
        const scalar betaOverbarf = betaOverbarIn[faceI];
        const scalar betaBlend = betaPlus*betaMinus/betaOverbarf;
        const scalar betaMinusOverbar = betaMinus/betaOverbarf;
        const scalar betaPlusOverbar = betaPlus/betaOverbarf;
        const scalar hJumpf = (2*wetOwner - 1)*hJumpIn[faceI];

        l[faceI] *= (1 - wetOwner)*betaMinusOverbar + wetOwner*betaPlusOverbar;
        u[faceI] *= wetOwner*betaMinusOverbar + (1 - wetOwner)*betaPlusOverbar;


        source[own] +=  betaBlend*cellVIn[own]*ownLs[faceI]*hJumpf;
        source[nei] -=  betaBlend*cellVIn[nei]*neiLs[faceI]*hJumpf;
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

            // Get references to interface related patch data
            const fvsPatchScalarField& pwetOwners =
                wetOwners.boundaryField()[patchI];
            const fvsPatchScalarField& pbetaOverbar =
                betaOverbar.boundaryField()[patchI];
            const fvsPatchScalarField& phJump = hJump.boundaryField()[patchI];
            const interfaceFvData::DynamicLabelList& pinterfaceFaces =
                boundaryInterfaceFaces[patchI];
            const vectorField& pownLs = ownLs.boundaryField()[patchI];

            // Get reference to face cells
            const unallocLabelList& pFaceCells =
                mesh.boundary()[patchI].faceCells();

            forAll(pinterfaceFaces, pifI)
            {
                // Get face and cell index
                const label& pfaceI = pinterfaceFaces[pifI];
                const label& pcellI = pFaceCells[pfaceI];

                // Get wet owner blending
                const scalar& pwetOwner = pwetOwners[pfaceI];

                // Helper variables
                const scalar pbetaCell = betaCell[pcellI];
                const scalar jumpCorr =
                betaPlus*betaMinus/pbetaOverbar[pfaceI]/pbetaCell;
                const scalar pbetaBlend =
                betaPlus*betaMinus/pbetaOverbar[pfaceI];

                // Remove ordinary 1/rho conditions to upper and lower coeffs and
                // add jump conditions
                pcL[pfaceI] *= jumpCorr; 
                pcU[pfaceI] *= jumpCorr;

                // Remove ordinary diagonal contribution and add jump condition
                d[pcellI] += cellVIn[pcellI]*pownLs[pfaceI]*pbetaCell;
                d[pcellI] -= cellVIn[pcellI]*pownLs[pfaceI]*pbetaBlend;

                // Add jump condition contribution to the source
                source[pcellI] += pbetaBlend*cellVIn[pcellI]*pownLs[pfaceI]*
                    (2*pwetOwner - 1)*phJump[pfaceI];
            }
        }

    }

    // Construct diagonals
    forAll(neighbour, faceI)
    {
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];
        // Caution - this is NOT negSumDiag(). VV, 17/July/2014 (IG, 6/Oct/2015)
        d[own] -= u[faceI];
        d[nei] -= l[faceI];
    }

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
