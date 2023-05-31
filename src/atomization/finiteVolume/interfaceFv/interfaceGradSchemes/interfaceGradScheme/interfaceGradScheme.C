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

#include "interfaceGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceGradScheme::~interfaceGradScheme()
{}


// * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * * * //

void interfaceGradScheme::interfaceCorrectBoundaryConditions
(
    const volScalarField& pd,
    volVectorField& gradpd,
    const interfaceFvData& intFvData
)
{
    // Get reference to interface finite volume data
    const surfaceScalarField& wetOwners = intFvData.wetOwnersOld();

    const scalar betaPlus = intFvData.betaPlus();
    const scalar betaMinus = intFvData.betaMinus();

    // Correct non coupled boundaries with snGrad of pd at the boundary
    forAll (pd.boundaryField(), patchI)
    {
        const fvPatchScalarField& ppd = pd.boundaryField()[patchI];

        if (!ppd.coupled())
        {
            const fvsPatchScalarField& pwetOwners =
                wetOwners.boundaryField()[patchI];

            const vectorField n = pd.mesh().boundary()[patchI].nf();
            const scalarField sngpd = ppd.snGrad();

            vectorField& pgpd = gradpd.boundaryField()[patchI];

            forAll(pgpd, pfaceI)
            {
                // Get normal vector
                const vector& nf = n[pfaceI];

                // Get gradient to correct
                vector& gradf = pgpd[pfaceI];

                // Bugfix: subtract the existing gradient in the normal
                // direction without 1/rho (betaPlus or betaMinus) since it's
                // already inside. VV, 2/Jul/2019
                gradf -= nf*(nf & gradf);

                // Now add the normal component from snGrad which needs to be
                // weighted with 1/rho
                gradf += nf*sngpd[pfaceI]*
                (
                    pwetOwners[pfaceI]*betaPlus
                  + (1.0 - pwetOwners[pfaceI]*betaMinus)
                );
            }
        }
    }

    // Evaluate coupled patches to correct patchNeighbourField
    gradpd.boundaryField().updateCoupledPatchFields();
}


tmp
<
    BlockLduSystem<vector, vector>
>
interfaceGradScheme::fvmGrad(const volScalarField& pd) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> interfaceGradScheme::fvmGrad\n"
        "(\n"
        " volScalarField& "
        ")\n"
    )   << "Implicit gradient operator currently defined only for Gauss linear "
        << "and leastSquares (cell and face limiters are optional)."
        << abort(FatalError);

    typedef typename outerProduct<vector, scalar>::type GradType;

    tmp<BlockLduSystem<vector, GradType> > tbs
    (
        new BlockLduSystem<vector, GradType>(pd.mesh())
    );

    return tbs;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
