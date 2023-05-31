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

Class
    Foam::fv::interfaceUncorrectedSnGrad

Description
    Simple interface corrected central-difference snGrad scheme without
    non-orthogonal correction.

SourceFiles
    interfaceUncorrectedSnGrad.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceUncorrectedSnGrad_H
#define interfaceUncorrectedSnGrad_H

#include "interfaceSnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class interfaceUncorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

class interfaceUncorrectedSnGrad
:
    public interfaceSnGradScheme
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const interfaceUncorrectedSnGrad&);


public:

    //- Runtime type information
    TypeName("interfaceUncorrected");


    // Constructors

        //- Construct from mesh
        interfaceUncorrectedSnGrad(const fvMesh& mesh)
        :
            interfaceSnGradScheme(mesh)
        {}


        //- Construct from mesh and data stream
        interfaceUncorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            interfaceSnGradScheme(mesh)
        {}


    // Destructor

        virtual ~interfaceUncorrectedSnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const volScalarField& pd
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the interfaceUncorrectedSnGrad
        //  for the given field
        virtual tmp<surfaceScalarField> correction
        (
            const volScalarField& pd
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
