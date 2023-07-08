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
    Simple central-difference snGrad scheme without non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "interfaceUncorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

defineTypeNameAndDebug(interfaceUncorrectedSnGrad, 0);
snGradScheme<scalar>::addMeshConstructorToTable<interfaceUncorrectedSnGrad>
    addInterfaceUncorrectedSnGradScalarMeshConstructorToTable_;


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interfaceUncorrectedSnGrad::~interfaceUncorrectedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<surfaceScalarField> interfaceUncorrectedSnGrad::correction
(
    const volScalarField& pd
) const
{
    notImplemented
    (
        "interfaceUncorrectedSnGrad::correction"
        "("
        "    volScalarField& pd"
        "    const interfaceFvData& intFvData"
        ")"
    );

    return tmp<surfaceScalarField>(nullptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
