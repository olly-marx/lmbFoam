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
    Foam::curvatureModel

Description
    Abstract base class that provides run-time selectable curvature calculation
    models.

SourceFiles
    curvatureModel.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef curvatureModel_H
#define curvatureModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class curvatureModel Declaration
\*---------------------------------------------------------------------------*/

class curvatureModel
{
    // Private Data Members

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Dictionary containing all the controls
        const dictionary coeffDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        curvatureModel(const curvatureModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const curvatureModel&) = delete;


protected:

    // Protected member functions

        // Access functions for derived classes

            //- Const access to the fvMesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Const access to the coefficient dictionary
            const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


public:

    //- Runtime type information
    TypeName("curvatureModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvatureModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        curvatureModel(const fvMesh& mesh, const dictionary& dict);


    // Selectors

        //- Return an autoPtr to the selected curvatureModel
        static autoPtr<curvatureModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~curvatureModel();


    // Member Functions

        //- Calculate and return mean curvature field given alpha field
        virtual tmp<volScalarField> kappa
        (
            const volScalarField& alpha
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
