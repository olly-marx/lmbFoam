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
    Foam::fv::interfaceGradScheme

Description
    Abstract base class for gradient schemes with interface correction. Intended
    only for the grad(pd)/rho term in the momentum equation for two phase flows
    (and for non-orthogonal correction for interfaceLaplacian and
    interfaceSnGrad). 1/rho jump is taken into account using interfaceFvData.
    interfaceGrad uses old interface data. User of the class has to update old
    interface data accordingly prior to the call.

SourceFiles
    interfaceGradScheme.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceGradScheme_H
#define interfaceGradScheme_H

#include "gradScheme.H"
#include "interfaceFvData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class interfaceGradScheme Declaration
\*---------------------------------------------------------------------------*/

class interfaceGradScheme
:
    public gradScheme<scalar>
{
    // Private Member Functions

        //- Disallow copy construct
        interfaceGradScheme(const interfaceGradScheme&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceGradScheme&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Constructors

        //- Construct from mesh
        interfaceGradScheme(const fvMesh& mesh)
        :
            gradScheme<scalar>(mesh)
        {}


    // Destructor

        virtual ~interfaceGradScheme();


    // Member Functions


        //- Calculate and return the grad of the given field
        virtual tmp<volVectorField> calcGrad
        (
            const volScalarField& pd,
            const word& name
        ) const = 0;


        //- Return the BlockLduSystem corresponding to the implicit grad
        // discretization. For block coupled systems.
        virtual tmp
        <
            BlockLduSystem<vector, vector> 
        > fvmGrad
        (
            const volScalarField& pd
        ) const;


        //- Correct the boundary values of the gradient using the patchField
        //  snGrad functions.
        static void interfaceCorrectBoundaryConditions
        (
            const volScalarField& pd,
            volVectorField& gradpd,
            const interfaceFvData& intFvData
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
