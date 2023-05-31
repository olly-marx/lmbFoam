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
    Foam::fv::interfaceGaussGrad

Description
    Basic second-order gradient scheme using face-interpolation and Gauss'
    theorem with second order correction at the interface. 1/rho jump
    taken into account using interfaceFvData.

SourceFiles
    interfaceGuassGrad.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceGaussGrad_H
#define interfaceGaussGrad_H

#include "interfaceGradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class interfaceGaussGrad Declaration
\*---------------------------------------------------------------------------*/

class interfaceGaussGrad
:
    public interfaceGradScheme
{
    // Private data

        //- Interpolation scheme to interpolate the pressure field onto faces.
        //  Interpolation is carried manually inside grad function for easier
        //  special treatment of interface faces
        tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceGaussGrad(const interfaceGaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceGaussGrad&);


public:

    //- Runtime type information
    TypeName("interfaceGauss");


    // Constructors

        //- Construct from mesh
        interfaceGaussGrad(const fvMesh& mesh)
        :
            interfaceGradScheme(mesh),
            tinterpScheme_(new linear<scalar>(mesh))
        {}

        //- Construct from Istream
        interfaceGaussGrad(const fvMesh& mesh, Istream& is)
        :
            interfaceGradScheme(mesh),
            tinterpScheme_(nullptr)
        {
            if (is.eof())
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<scalar> >
                    (
                        new linear<scalar>(mesh)
                    );
            }
            else
            {
                tinterpScheme_ =
                    tmp<surfaceInterpolationScheme<scalar> >
                    (
                        surfaceInterpolationScheme<scalar>::New(mesh, is)
                    );
            }
        }


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
