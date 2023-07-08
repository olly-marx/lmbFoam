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
    Foam::fv::interfaceSnGradScheme

Description
    Abstract base class for snGrad schemes with interface correction. Intended
    only for the grad(pd)/rho term in the momentum equation for two phase
    flows.
    1/rho jump is taken into account using interfaceFvData.
    interfaceSnGrad uses old interface data. User of the class has to
    update old interface data accordingly prior to the call.

SourceFiles
    interfaceSnGradScheme.C

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceSnGradScheme_H
#define interfaceSnGradScheme_H

#include "snGradScheme.H"
#include "interfaceFvData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class interfaceSnGradScheme Declaration
\*---------------------------------------------------------------------------*/

class interfaceSnGradScheme
:
    public snGradScheme<scalar>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const interfaceSnGradScheme&);


public:

    // Constructors

        //- Construct from mesh
        interfaceSnGradScheme(const fvMesh& mesh)
        :
            snGradScheme<scalar>(mesh)
        {}


    // Destructor

        virtual ~interfaceSnGradScheme();


    // Member Functions

        //- Return the snGrad of the given cell field with the given
        //  deltaCoeffs
        static tmp<surfaceScalarField> interfaceSnGrad
        (
            const volScalarField& pd,
            const tmp<surfaceScalarField>& tdeltaCoeffs,
            const interfaceFvData& intFvData,
            const word& snGradName = "snGrad"
        );

        //- Correct orthogonal part using jump conditions
        static void orthogonalInterfaceCorrection
        (
            surfaceScalarField& orthogonalSnGrad,
            const tmp<surfaceScalarField>& tdeltaCoeffs,
            const interfaceFvData& intFvData
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const volScalarField& pd
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the snGrad
        //  for the given field
        virtual tmp<surfaceScalarField> correction
        (
            const volScalarField& pd
        ) const
        {
            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the snGrad of the given cell field
        //  with explicit correction
        virtual tmp<surfaceScalarField> snGrad
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
