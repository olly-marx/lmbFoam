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
    Foam::interpolationSmoothedLevelSet

Description
    Calculates the curvature based on the level set field. The level set field
    (signed distance functions) is calculated using interfaceFacesWave class
    that calculates the distance from the interface faces using MeshWave
    procedure. The level set field is optionally smoothed using volume-to-point
    and then point-to-volume interpolation.

SourceFiles
    interpolationSmoothedLevelSet.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interpolationSmoothedLevelSet_H
#define interpolationSmoothedLevelSet_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class interpolationSmoothedLevelSet Declaration
\*---------------------------------------------------------------------------*/

class interpolationSmoothedLevelSet
:
    public curvatureModel
{
    // Private data

        //- Number of level set smoothing iterations
        label nSmoothIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interpolationSmoothedLevelSet
        (
            const interpolationSmoothedLevelSet&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const interpolationSmoothedLevelSet&) = delete;


public:

    //- Runtime type information
    TypeName("interpolationSmoothedLevelSet");


    // Constructors

        //- Construct from components
        interpolationSmoothedLevelSet
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interpolationSmoothedLevelSet();


    // Member functions

        //- Calculate and return mean curvature field.
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
