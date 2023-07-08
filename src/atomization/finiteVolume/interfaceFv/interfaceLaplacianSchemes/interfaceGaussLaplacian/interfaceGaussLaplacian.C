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

#include "interfaceGaussLaplacian.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "interfaceSnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceGaussLaplacian, 0);
laplacianScheme<scalar, scalar>::addIstreamConstructorToTable
<
    interfaceGaussLaplacian
>
addInterfaceGaussLaplacianScalarIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void interfaceGaussLaplacian::checkSnGradScheme() const
{
    if (!isA<interfaceSnGradScheme>(this->tsnGradScheme_()))
    {
        FatalErrorIn
        (
            "interfaceGaussLaplacian::checkSnGradScheme()"
        ) << "snGrad scheme for pressure laplacian is not"
          << " interface corrected." << nl
          << "Please use interface corrected snGrad scheme for pressure"
          << "laplacian in order to use the Ghost Fluid Method." << endl
          << exit(FatalError);
    }
}


tmp<fvScalarMatrix>
interfaceGaussLaplacian::fvmInterfaceLaplacianUncorrected
(
    const surfaceScalarField& rAUfMagSf,
    const volScalarField& pd,
    const interfaceFvData& intFvData
)
{
    tmp<surfaceScalarField> tdeltaCoeffs =
        this->tsnGradScheme_().deltaCoeffs(pd);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            pd,
            deltaCoeffs.dimensions()*rAUfMagSf.dimensions()
           *pd.dimensions()/intFvData.rhoPlus().dimensions()
        )
    );
    fvScalarMatrix& fvm = tfvm();

    fvm.upper() = deltaCoeffs.internalField()*rAUfMagSf.internalField();
    // negSumDiag will be called after the interface jump conditions

    forAll(fvm.psi().boundaryField(), patchI)
    {
        const fvPatchScalarField& psf = fvm.psi().boundaryField()[patchI];
        const fvsPatchScalarField& patchrAUf =
            rAUfMagSf.boundaryField()[patchI];

        fvm.internalCoeffs()[patchI] = patchrAUf*psf.gradientInternalCoeffs();
        fvm.boundaryCoeffs()[patchI] = -patchrAUf*psf.gradientBoundaryCoeffs();
    }

    return tfvm;
}


void interfaceGaussLaplacian::correctPressureEquationAndFlux
(
    fvScalarMatrix& pdEqn,
    surfaceScalarField& phiJump,
    const interfaceFvData& intFvData
)
{
    // Get reference to const mesh data
    const fvMesh& mesh = this->mesh();
    const unallocLabelList& neighbour = mesh.neighbour();

    // STAGE 1: density (1/rho) jump condition

    // Only upper coefficients needed because the resulting matrix will be
    // symmetric and the diagonal can be reconstructed using negSumDiag(). The
    // source is at this stage zero (no non-orthogonal correction is added yet).
    scalarField& upper = pdEqn.upper();

    // Add 1/rho in upper coefficients, taking the jump in rho into account
    // Get necessary constants from interfaceFvData object
    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Get reference to wet owners (wet owner = 1, dry owner = 0)
    const surfaceScalarField& wetOwners = intFvData.wetOwners();
    const scalarField& wetOwnersIn = wetOwners.internalField();

    // Scale upper matrix coefficients
    forAll(neighbour, faceI)
    {
        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Scale upper coefficients at this face
        upper[faceI] *= wetOwner*betaPlus + (1 - wetOwner)*betaMinus;
    }

    // Boundary contributions
    forAll(pdEqn.boundaryCoeffs(), patchI)
    {
        // Get wetOwners patch field
        const fvsPatchScalarField& pwetOwners =
            wetOwners.boundaryField()[patchI];

        const scalarField prhoBlend
        (
            pwetOwners*betaPlus + (1 - pwetOwners)*betaMinus
        );

        pdEqn.internalCoeffs()[patchI] *= prhoBlend;
        pdEqn.boundaryCoeffs()[patchI] *= prhoBlend;
    }

    // STAGE 2: pressure gradient (grad(pd)) jump condition

    // Get reference to internal flux field
    scalarField& phiJumpIn = phiJump.internalField();

    // Get references to needed interface finite volume data
    const surfaceScalarField& betaOverbar = intFvData.betaOverbar();
    const scalarField& betaOverbarIn = betaOverbar.internalField();
    const surfaceScalarField& hJump = intFvData.hydrostaticJump();
    const scalarField& hJumpIn = hJump.internalField();

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFaces();

    // Add jump conditions to internal interface faces
    forAll(interfaceFaces, ifI)
    {
        // Get face index from list of interface faces
        const label& faceI = interfaceFaces[ifI];

        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Helper variables
        const scalar& betaOverbarf = betaOverbarIn[faceI];
        const scalar betaMinusOverbar = betaMinus/betaOverbarf;
        const scalar betaPlusOverbar = betaPlus/betaOverbarf;

        // Get reference to off diagonal coefficient
        scalar& aUpper = upper[faceI];

        // Add pressure jump condition contributions
        // Source contribution is asymmetric - add it to the flux before final
        // tempering with upper coefficient. Note: negative sign (assumed
        // already on the RHS of equation)
        phiJumpIn[faceI] -= aUpper*hJumpIn[faceI]*
            (wetOwner*betaMinusOverbar - (1 - wetOwner)*betaPlusOverbar);

        // Off diagonal contribution - only upper because after this jump, the
        // matrix remains symmetric
        aUpper *= wetOwner*betaMinusOverbar + (1 - wetOwner)*betaPlusOverbar;
    }

    // Calculate diagonal using negSumDiag. No need to reset diagonal because it
    // is not calculated in fvmInterfaceLaplacianUncorrected() member
    pdEqn.negSumDiag();

    // Get reference to boundary interface faces
    const interfaceFvData::DynamicLabelListList& boundaryInterfaceFaces =
        intFvData.boundaryInterfaceFaces();

    // Add jump conditions to boundary interface faces
    forAll(boundaryInterfaceFaces, patchI)
    {
        const fvPatchScalarField& ppd = pdEqn.psi().boundaryField()[patchI];

        if (ppd.coupled())
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
            scalarField& pinternalCoeffs = pdEqn.internalCoeffs()[patchI];
            scalarField& pboundaryCoeffs = pdEqn.boundaryCoeffs()[patchI];
            scalarField& pphiJump = phiJump.boundaryField()[patchI];

            forAll(pinterfaceFaces, pifI)
            {
                // Get patch face index
                const label& pfaceI = pinterfaceFaces[pifI];

                // Get wet owner blending
                const scalar& pwetOwner = pwetOwners[pfaceI];

                // Helper variables
                const scalar& pbetaOverbarf = pbetaOverbar[pfaceI];
                const scalar pbetaMinusOverbar = betaMinus/pbetaOverbarf;
                const scalar pbetaPlusOverbar = betaPlus/pbetaOverbarf;

                // Add additional source term. Note: opposite sign because of
                // boundary coeffs.
                pphiJump[pfaceI] += pboundaryCoeffs[pfaceI]*phJump[pfaceI]*
                (
                    pwetOwner*pbetaMinusOverbar
                  - (1 - pwetOwner)*pbetaPlusOverbar
                );

                const scalar pbetaBlend = pwetOwner*pbetaMinusOverbar
                    + (1 - pwetOwner)*pbetaPlusOverbar;

                // Scale diagonal and source contributions.
                pboundaryCoeffs[pfaceI] *= pbetaBlend;
                pinternalCoeffs[pfaceI] *= pbetaBlend;
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015.
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<fvScalarMatrix>
interfaceGaussLaplacian::fvmLaplacian
(
    const surfaceScalarField& rAUf,
    const volScalarField& pd
)
{
    const fvMesh& mesh = this->mesh();

    // Calculate Sf/a
    const surfaceScalarField rAUfMagSf = rAUf*mesh.magSf();

    // Lookup interface data
    const interfaceFvData& intFvData =
        mesh.lookupObject<interfaceFvData>("interfaceFvData");

    // Refactorization: phiJump iis now lumped into faceFluxCorrection

    // Create face flux correction for jump conditions and non-orthogonal
    // corrected schemes
    tmp<surfaceScalarField> tfaceFluxCorrection
    (
        new surfaceScalarField
        (
            IOobject
            (
                "faceFluxCorrection(" + pd.name() + ')',
                pd.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                mesh.deltaCoeffs().dimensions()*rAUfMagSf.dimensions()
               *pd.dimensions()/intFvData.rhoPlus().dimensions(),
                0
            )
        )
    );
    surfaceScalarField& fluxCorrection = tfaceFluxCorrection();

    // Calculate usual part of the matrix without 1/rho and pd jump
    tmp<fvScalarMatrix> tfvm = fvmInterfaceLaplacianUncorrected
    (
        rAUfMagSf,
        pd,
        intFvData
    );
    fvScalarMatrix& fvm = tfvm();

    // Correct the pressure equation and the flux due to interface jump
    // conditions (1/rho and grad(pd))
    correctPressureEquationAndFlux(fvm, fluxCorrection, intFvData);

    // Non-orthogonal correction from interfaceSnGradScheme
    if (this->tsnGradScheme_().corrected())
    {
        fluxCorrection += rAUfMagSf*this->tsnGradScheme_().correction(pd);
    }

    // Add the jump condition flux and non-orthogonal correction to the matrix.
    // Note: fvMatrix::operator+= subtracts from source. Jump condition flux
    // sign is now changed
    fvm += fvc::div(fluxCorrection);

    // Add the non-orthogonal correction in the face flux correction pointer
    if (mesh.schemesDict().fluxRequired(pd.name()))
    {
        fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
    }

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
