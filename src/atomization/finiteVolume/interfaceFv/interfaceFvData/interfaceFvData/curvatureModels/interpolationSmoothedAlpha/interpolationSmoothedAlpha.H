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
    Foam::interpolationSmoothedAlpha

Description
    Calculates the curvature based on gradient of alpha field. The gradient is
    smoothed by interpolating from cell centres to points and then back from
    points to cell centres a number of times.

SourceFiles
    interpolationSmoothedAlpha.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interpolationSmoothedAlpha_H
#define interpolationSmoothedAlpha_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class interpolationSmoothedAlpha Declaration
\*---------------------------------------------------------------------------*/

class interpolationSmoothedAlpha
:
    public curvatureModel
{
    // Private data

        //- Number of smoothing iterations
        label nSmoothIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interpolationSmoothedAlpha
        (
            const interpolationSmoothedAlpha&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=
        (
            const interpolationSmoothedAlpha&
        ) = delete;


public:

    //- Runtime type information
    TypeName("interpolationSmoothedAlpha");


    // Constructors

        //- Construct from components
        interpolationSmoothedAlpha
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interpolationSmoothedAlpha();


    // Member functions

        //- Calculate and return mean curvature field given alpha field
        virtual tmp<volScalarField> kappa
        (
            const volScalarField& alpha
        ) const final;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
