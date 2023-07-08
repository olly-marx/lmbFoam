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
    Foam::fv::interfaceLimitedSnGrad

Description
    Central-difference interface corrected snGrad scheme with limited
    non-orthogonal correction.

    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    interfaceUncorrectedSnGrad, when set to 1 the full correction is applied
    and the scheme behaves as interfaceCorrectedSnGrad and when set to 0.5 the
    limiter is calculated such that the non-orthogonal contribution does not
    exceed the orthogonal part.

SourceFiles
    interfaceLimitedSnGrad.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceLimitedSnGrad_H
#define interfaceLimitedSnGrad_H

#include "interfaceSnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class interfaceLimitedSnGrad Declaration
\*---------------------------------------------------------------------------*/

class interfaceLimitedSnGrad
:
    public interfaceSnGradScheme
{
    // Private data

        scalar limitCoeff_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const interfaceLimitedSnGrad&);


public:

    //- Runtime type information
    TypeName("interfaceLimited");


    // Constructors

        //- Construct from mesh
        interfaceLimitedSnGrad(const fvMesh& mesh)
        :
            interfaceSnGradScheme(mesh)
        {}


        //- Construct from mesh and data stream
        interfaceLimitedSnGrad(const fvMesh& mesh, Istream& is)
        :
            interfaceSnGradScheme(mesh),
            limitCoeff_(readScalar(is))
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorIn
                (
                    "interfaceLimitedSnGrad(const fvMesh& mesh, Istream& is): ",
                    is
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }


    // Destructor

        virtual ~interfaceLimitedSnGrad();


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
            return !this->mesh().orthogonal();
        }

        //- Return the explicit correction to the interfaceLimitedSnGrad
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
