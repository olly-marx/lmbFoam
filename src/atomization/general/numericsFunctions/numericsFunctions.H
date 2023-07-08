/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

InNamespace
    Foam::numericsFunctions

Description
    General useful functions related to different aspects of numerics in
    atomization simulations

SourceFiles
    numericsFunctions.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef numericsFunctions_H
#define numericsFunctions_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace numericsFunctions
{

//- Limit the internal vector field given field to limit and dictionary to read
//  the limit from
//  Parameters:
//  1. Velocity field to limit
//  2. Dictionary to read the limitMag value from
void limitVelocityField
(
    volVectorField& UField,
    const dictionary& dict
);


//- Extrapolate velocity field from water side near the interface. A similar
//  treatment to this one is used in CFDShip-Iowa in their single phase solver
void extrapolateInterfaceVelocity
(
    volVectorField& UField,
    const dictionary& dict
);


} // End namespace numericsFunction

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
