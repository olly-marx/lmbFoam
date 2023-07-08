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
    Foam::fv::interfaceLeastSquaresGrad

Description
    Second-order gradient scheme using least-squares with correction at the
    interface. 1/rho jump taken into account using interfaceFvData.

SourceFiles
    interfaceLeastSquaresGrad.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceLeastSquaresGrad_H
#define interfaceLeastSquaresGrad_H

#include "interfaceGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class interfaceLeastSquaresGrad Declaration
\*---------------------------------------------------------------------------*/

class interfaceLeastSquaresGrad
:
    public interfaceGradScheme
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceLeastSquaresGrad(const interfaceLeastSquaresGrad&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceLeastSquaresGrad&);


public:

    //- Runtime type information
    TypeName("interfaceLeastSquares");


    // Constructors

        //- Construct from mesh
        interfaceLeastSquaresGrad(const fvMesh& mesh)
        :
            interfaceGradScheme(mesh)
        {}

        //- Construct from Istream
        interfaceLeastSquaresGrad(const fvMesh& mesh, Istream&)
        :
            interfaceGradScheme(mesh)
        {}


    // Member Functions

        //- Calculate and return grad(pd)
        virtual tmp<volVectorField> calcGrad
        (
            const volScalarField& pd,
            const word& name
        ) const;

        //- Return the BlockLduSystem corresponding to the implicit grad
        // discretization. For block coupled systems.
        virtual tmp
        <
            BlockLduSystem<vector, vector>
        > fvmGrad
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
