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

#include "interfaceCellLimitedGrad.H"
#include "interfaceGaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceCellLimitedGrad, 0);
gradScheme<scalar>::addIstreamConstructorToTable<interfaceCellLimitedGrad>
    addInterfaceCellLimitedScalarIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void interfaceCellLimitedGrad::limitFace
(
    scalar& limiter,
    const scalar& maxDelta,
    const scalar& minDelta,
    const scalar& extrapolate
)
{
    if (extrapolate > maxDelta + VSMALL)
    {
        limiter = min(limiter, maxDelta/extrapolate);
    }
    else if (extrapolate < minDelta - VSMALL)
    {
        limiter = min(limiter, minDelta/extrapolate);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volVectorField> interfaceCellLimitedGrad::calcGrad
(
    const volScalarField& pd,
    const word& name
) const
{
    const fvMesh& mesh = pd.mesh();

    tmp<volVectorField> tGrad = basicInterfaceGradScheme_().calcGrad(pd, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    volVectorField& g = tGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalarField maxpd(pd.internalField());
    scalarField minpd(pd.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar pdOwn = pd[own];
        scalar pdNei = pd[nei];

        maxpd[own] = max(maxpd[own], pdNei);
        minpd[own] = min(minpd[own], pdNei);

        maxpd[nei] = max(maxpd[nei], pdOwn);
        minpd[nei] = min(minpd[nei], pdOwn);
    }

    const volScalarField::GeometricBoundaryField& bpd = pd.boundaryField();

    forAll(bpd, patchi)
    {
        const fvPatchScalarField& ppd = bpd[patchi];

        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();

        if (ppd.coupled())
        {
            scalarField ppdNei = ppd.patchNeighbourField();

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar pdNei = ppdNei[pFacei];

                maxpd[own] = max(maxpd[own], pdNei);
                minpd[own] = min(minpd[own], pdNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar pdNei = ppd[pFacei];

                maxpd[own] = max(maxpd[own], pdNei);
                minpd[own] = min(minpd[own], pdNei);
            }
        }
    }

    maxpd -= pd;
    minpd -= pd;

    if (k_ < 1.0)
    {
        scalarField maxMinpd = (1.0/k_ - 1.0)*(maxpd - minpd);
        maxpd += maxMinpd;
        minpd -= maxMinpd;

        //maxpd *= 1.0/k_;
        //minpd *= 1.0/k_;
    }


    // create limiter
    scalarField limiter(pd.internalField().size(), 1.0);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxpd[own],
            minpd[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxpd[nei],
            minpd[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bpd, patchi)
    {
        const unallocLabelList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxpd[own],
                minpd[own],
                (pCf[pFacei] - C[own]) & g[own]
            );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << pd.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    g.internalField() *= limiter;

    // Correct boundary conditions
    g.correctBoundaryConditions();
    interfaceGradScheme::interfaceCorrectBoundaryConditions
    (
        pd,
        g,
        mesh.lookupObject<interfaceFvData>("interfaceFvData")
    );

    return tGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
