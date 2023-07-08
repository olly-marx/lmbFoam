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

#include "interfaceSnGradScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceSnGradScheme::~interfaceSnGradScheme()
{}


// * * * * * * * * * * * Static member Functions * * * * * * * * * * * * * * //

tmp<surfaceScalarField> interfaceSnGradScheme::interfaceSnGrad
(
    const volScalarField& pd,
    const tmp<surfaceScalarField>& tdeltaCoeffs,
    const interfaceFvData& intFvData,
    const word& snGradName
)
{
    // Get reference to mesh
    const fvMesh& mesh = pd.mesh();

    // Construct surfaceScalarField
    tmp<surfaceScalarField> tssf
    (
        new surfaceScalarField
        (
            IOobject
            (
                snGradName + "(" + pd.name() + ')',
                pd.instance(),
                pd.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            pd.dimensions()
           *tdeltaCoeffs().dimensions()/intFvData.rhoPlus().dimensions()
        )
    );
    surfaceScalarField& ssf = tssf();

    // Get reference to difference factors
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

    // Get owner/neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get reference to internal pd field
    const scalarField& pdIn = pd.internalField();

    // Get reference to internal snGrad field
    scalarField& ssfIn = ssf.internalField();

    // Ordinary snGrad for internal faces
    forAll(neighbour, faceI)
    {
        ssfIn[faceI] =
            deltaCoeffs[faceI]*(pdIn[neighbour[faceI]] - pdIn[owner[faceI]]);
    }

    // Ordinary snGrad for boundary faces
    forAll(pd.boundaryField(), patchI)
    {
        ssf.boundaryField()[patchI] = pd.boundaryField()[patchI].snGrad();
    }

    // Correct orthogonal part using jump conditions
    orthogonalInterfaceCorrection(ssf, tdeltaCoeffs, intFvData);

    return tssf;
}


void interfaceSnGradScheme::orthogonalInterfaceCorrection
(
    surfaceScalarField& orthogonalSnGrad,
    const tmp<surfaceScalarField>& tdeltaCoeffs,
    const interfaceFvData& intFvData
)
{
    // Get reference to const mesh data
    const fvMesh& mesh = orthogonalSnGrad.mesh();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get reference to orthogonalSnGrad internalField
    scalarField& orthogonalSnGradIn = orthogonalSnGrad.internalField();

    // STAGE 1: density (1/rho) jump condition

    // Scale the orthogonal snGrad with 1/rho
    // Get necessary constants from interfaceFvData object
    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Get reference to wet owners (wet owner = 1, dry owner = 0)
    const surfaceScalarField& wetOwners = intFvData.wetOwnersOld();
    const scalarField& wetOwnersIn = wetOwners.internalField();

    // Internal face contributions
    forAll(neighbour, faceI)
    {
        // Get wet owner blending
        const scalar& wetOwner = wetOwnersIn[faceI];

        // Scale orthogonal snGrad at this face
        orthogonalSnGradIn[faceI] *=
            wetOwner*betaPlus + (1 - wetOwner)*betaMinus;
    }

    // Boundary contributions
    forAll(orthogonalSnGrad.boundaryField(), patchI)
    {
        // Get wetOwners patch field
        const fvsPatchScalarField& pwetOwners =
            wetOwners.boundaryField()[patchI];

        orthogonalSnGrad.boundaryField()[patchI] *=
            pwetOwners*betaPlus + (1 - pwetOwners)*betaMinus;
    }

    // STAGE 2: pressure gradient jump condition

    // Get references to needed interface finite volume data
    const surfaceScalarField& betaOverbar = intFvData.betaOverbarOld();
    const scalarField& betaOverbarIn = betaOverbar.internalField();
    const surfaceScalarField& hJump = intFvData.hydrostaticJumpOld();
    const scalarField& hJumpIn = hJump.internalField();

    // Get reference to the list of interfaceFaces
    const interfaceFvData::DynamicLabelList& interfaceFaces =
        intFvData.interfaceFacesOld();

    // Get reference to internal deltaCoeffs field for additional hydrostatic
    // jump
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

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

        // First scale the orthogonalSnGrad field...
        orthogonalSnGradIn[faceI] *=
            wetOwner*betaMinusOverbar + (1 - wetOwner)*betaPlusOverbar;

        // ... then add additional source due to jump condition
        orthogonalSnGradIn[faceI] += deltaCoeffs[faceI]*hJumpIn[faceI]*
        // Left for clarity even though:
        // betaMinusOverbar*betaPlus == betaPlusOverbar*betaMinus
        (
          - wetOwner*betaMinusOverbar*betaPlus
          + (1 - wetOwner)*betaPlusOverbar*betaMinus
        );
    }

    // Get reference to boundary interface faces
    const interfaceFvData::DynamicLabelListList& boundaryInterfaceFaces =
        intFvData.boundaryInterfaceFacesOld();

    // Add jump conditions to boundary interface faces
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
            const scalarField& pdeltaCoeffs = patch.deltaCoeffs();
            fvsPatchScalarField& porthSnGrad =
                orthogonalSnGrad.boundaryField()[patchI];

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

                // First scale the snGrad at the patch face...
                porthSnGrad[pfaceI] *= pwetOwner*pbetaMinusOverbar
                  + (1 - pwetOwner)*pbetaPlusOverbar;

                // ... then add additional source due to jump condition
                porthSnGrad[pfaceI] += pdeltaCoeffs[pfaceI]*phJump[pfaceI]*
                // Left for clarity even though:
                // pbetaMinusOverbar*betaPlus == pbetaPlusOverbar*betaMinus
                (
                  - pwetOwner*pbetaMinusOverbar*betaPlus
                  + (1 - pwetOwner)*pbetaPlusOverbar*betaMinus
                );
            }
        }

        // Do not consider interface faces on non coupled patches! Both
        // fixedValue or fixedGradient pd boundaries do not make sense for
        // patches that might have interface faces and consequently pressure
        // jumps at those faces. VV, 6/March/2015.
    }
}


// * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> interfaceSnGradScheme::snGrad
(
    const volScalarField& pd
) const
{
    // Get reference to interfaceFvData object
    const interfaceFvData& intFvData =
        mesh().lookupObject<interfaceFvData>("interfaceFvData");

    // Calculate interface corrected orthogonal part of the gradient
    tmp<surfaceScalarField> tsf = interfaceSnGrad
    (
        pd,
        deltaCoeffs(pd),
        intFvData
    );

    // Calculate interface corrected non-orthogonal correction part of the
    // gradient
    if (corrected())
    {
        tsf() += correction(pd);
    }

    return tsf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
