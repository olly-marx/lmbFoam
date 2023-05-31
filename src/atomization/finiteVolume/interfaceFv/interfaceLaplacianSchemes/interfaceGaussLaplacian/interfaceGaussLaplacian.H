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
    Foam::fv::interfaceGaussLaplacian

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem with
    correction at the interface. Intended only for the implicit discretisation
    of the pressure equation in two phase flows. 1/rho taken into account from
    interfaceFvData class.

SourceFiles
    interfaceGaussLaplacian.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceGaussLaplacian_H
#define interfaceGaussLaplacian_H

#include "laplacianScheme.H"
#include "interfaceFvData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class interfaceGaussLaplacian Declaration
\*---------------------------------------------------------------------------*/

class interfaceGaussLaplacian
:
    public laplacianScheme<scalar, scalar>
{
    // Private Member Functions

        // Check that the selected snGrad scheme is interface corrected
        void checkSnGradScheme() const;

        //- Calculate uncorrected part of the pressure laplacian
        tmp<fvScalarMatrix> fvmInterfaceLaplacianUncorrected
        (
            const surfaceScalarField& rAUfMagSf,
            const volScalarField& pd,
            const interfaceFvData& intFvData
        );

        //- Embed jump conditions into the pressure equation and flux
        void correctPressureEquationAndFlux
        (
            fvScalarMatrix& pdEqn,
            surfaceScalarField& phi,
            const interfaceFvData& intFvData
        );

        //- Disallow default bitwise copy construct
        interfaceGaussLaplacian(const interfaceGaussLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceGaussLaplacian&);


public:

    //- Runtime type information
    TypeName("interfaceGauss");


    // Constructors

        //- Construct from mesh
        interfaceGaussLaplacian(const fvMesh& mesh)
        :
            laplacianScheme<scalar, scalar>(mesh)
        {
            checkSnGradScheme();
        }

        //- Construct from Istream
        interfaceGaussLaplacian(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<scalar, scalar>(mesh, is)
        {
            checkSnGradScheme();
        }

        //- Construct from mesh, interpolation and snGradScheme schemes
        interfaceGaussLaplacian
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<scalar> >& igs,
            const tmp<snGradScheme<scalar> >& sngs
        )
        :
            laplacianScheme<scalar, scalar>(mesh, igs, sngs)
        {
            checkSnGradScheme();
        }


    // Destructor

        virtual ~interfaceGaussLaplacian()
        {}


    // Member Functions

        // Implicit laplacian scheme for pressure equation only
        tmp<fvScalarMatrix> fvmLaplacian
        (
            const surfaceScalarField& rAUf,
            const volScalarField& pd
        );

        // Explicit laplacian schemes not implemented at the moment
        tmp<volScalarField> fvcLaplacian
        (
            const surfaceScalarField& rAUf,
            const volScalarField& pd
        )
        {
            notImplemented
            (
                "interfaceGaussLaplacian::fvcLaplacian"
            );

            // Dummy return to prevent compiler warnings. VV, 1/Apr/2015.
            return tmp<volScalarField>(nullptr);
        }

        tmp<volScalarField> fvcLaplacian
        (
            const volScalarField& pd
        )
        {
            notImplemented
            (
                "interfaceGaussLaplacian::fvcLaplacian"
            );

            // Dummy return to prevent compiler warnings. VV, 1/Apr/2015.
            return tmp<volScalarField>(nullptr);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
