/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Class
    Foam::sphereDistance

Description
    Distance class to the free surface defined by a sphere

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef sphereDistance_H
#define sphereDistance_H

#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sphereDistance Declaration
\*---------------------------------------------------------------------------*/

class sphereDistance
{
    // Private data

        //- Sphere centre
        const point centre_;

        //- Sphere radius
        const scalar radius_;


public:

    // Constructors

        //- Construct from components
        inline sphereDistance
        (
            const point& centre,
            const scalar radius
        )
        :
            centre_(centre),
            radius_(radius)
        {}


    //- Destructor
    ~sphereDistance() = default;


    // Member Functions

        // Iterate distance
        inline bool iterateDistance() const
        {
            return false;
        }

        // Return distance to the point
        inline scalar distance(const point& p) const
        {
            return mag(p - centre_) - radius_;
        }

        // Return distance to the point field
        inline tmp<scalarField> distance(const pointField& p) const
        {
            return mag(p - centre_) - radius_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
