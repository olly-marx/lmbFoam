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
    Foam::fv::interfaceFaceLimitedGrad

Description
    interfaceFaceLimitedGrad gradient scheme applied to a runTime selected base
    gradient scheme with interface corrections.

    The scalar limiter based on limiting the extrapolated face values
    between the face-neighbour cell values and is applied to all components
    of the gradient.

    Jump conditions are taken into account in interface grad schemes.
    Hence, the limiter is safe to use in its original form applied to
    interfaceGradSchemes.
    gradScheme::correctBoundaryCondition() at the end of grad
    calculation must not be called.

SourceFiles
    interfaceFaceLimitedGrad.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceFaceLimitedGrad_H
#define interfaceFaceLimitedGrad_H

#include "interfaceGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class interfaceFaceLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

class interfaceFaceLimitedGrad
:
    public interfaceGradScheme
{
    // Private Data

        //- Basic (unlimited) interface corrected grad scheme
        tmp<gradScheme<scalar> > basicInterfaceGradScheme_;

        //- Limiter coefficient
        const scalar k_;


    // Private Member Functions

        inline void limitFace
        (
            scalar& limiter,
            const scalar& maxDelta,
            const scalar& minDelta,
            const scalar& extrapolate
        ) const;


        //- Disallow default bitwise copy construct
        interfaceFaceLimitedGrad(const interfaceFaceLimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceFaceLimitedGrad&);


public:

    //- RunTime type information
    TypeName("interfaceFaceLimited");


    // Constructors

        //- Construct from mesh and schemeData
        interfaceFaceLimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            interfaceGradScheme(mesh),
            basicInterfaceGradScheme_
            (
                gradScheme<scalar>::New(mesh, schemeData)
            ),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorIn
                (
                    "interfaceFaceLimitedGrad"
                    "(const fvMesh& mesh, Istream& schemeData)",
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }

            // Sanity check the type
            if (!isA<interfaceGradScheme>(basicInterfaceGradScheme_()))
            {
                FatalIOErrorIn
                (
                    "interfaceFaceLimitedGrad"
                    "(const fvMesh& mesh, Istream& schemeData)",
                    schemeData
                )   << "Specified underlying basic gradient scheme: "
                    << basicInterfaceGradScheme_().type()
                    << nl
                    << "Make sure to use one of the interface corrected"
                    << " gradient schemes as basic grad scheme."
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        // - Calculate and return grad(pd)
        virtual tmp<volVectorField> calcGrad
        (
            const volScalarField& pd,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
