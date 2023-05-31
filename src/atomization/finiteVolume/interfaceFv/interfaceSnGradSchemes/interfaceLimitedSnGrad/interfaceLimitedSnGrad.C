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

Description
    Simple central-difference snGrad scheme with limited non-orthogonal
    correction.

\*---------------------------------------------------------------------------*/

#include "interfaceLimitedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interfaceCorrectedSnGrad.H"
#include "interfaceUncorrectedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceLimitedSnGrad, 0);
snGradScheme<scalar>::addMeshConstructorToTable<interfaceLimitedSnGrad>
    addInterfaceLimitedSnGradScalarMeshConstructorToTable_;


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceLimitedSnGrad::~interfaceLimitedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> interfaceLimitedSnGrad::correction
(
    const volScalarField& pd
) const
{
    // Get correction field from corrected scheme
    surfaceScalarField corr
    (
        interfaceCorrectedSnGrad(this->mesh()).correction(pd)
    );

    // Get reference to interfaceFvData object
    const interfaceFvData& intFvData =
        this->mesh().lookupObject<interfaceFvData>("interfaceFvData");

    // Calculate limiter from fully corrected and uncorrected parts
    surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag
            (
                interfaceSnGradScheme::interfaceSnGrad
                (
                    pd,
                    deltaCoeffs(pd),
                    intFvData,
                    "orthSnGrad"
                )
            )/
            (
                (1 - limitCoeff_)*mag(corr)
              + dimensionedScalar("small", corr.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (fv::debug)
    {
        Info<< "interfaceLimitedSnGrad :: limiter min: "
            << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
    }

    return limiter*corr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
