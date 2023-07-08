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
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "interfaceCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "gradScheme.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceCorrectedSnGrad, 0);
snGradScheme<scalar>::addMeshConstructorToTable<interfaceCorrectedSnGrad>
    addInterfaceCorrectedSnGradScalarMeshConstructorToTable_;


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceCorrectedSnGrad::~interfaceCorrectedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> interfaceCorrectedSnGrad::correction
(
    const volScalarField& pd
) const
{
    // Get reference to mesh
    const fvMesh& mesh = this->mesh();

    // Construct surface scalar field containing non-orthogonal correction
    tmp<surfaceScalarField> tssf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "snGradCorr(" + pd.name() + ')',
                pd.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            pd.dimensions()
           *mesh.deltaCoeffs().dimensions()/dimDensity
        )
    );
    surfaceScalarField& ssf = tssf();

    ssf =
        mesh.correctionVectors() & fvc::interpolate(fvc::grad(pd, ssf.name()));

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
