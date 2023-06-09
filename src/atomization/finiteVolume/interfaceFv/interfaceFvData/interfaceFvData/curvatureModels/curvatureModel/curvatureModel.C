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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "curvatureModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(curvatureModel, 0);
defineRunTimeSelectionTable(curvatureModel, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModel::curvatureModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    coeffDict_(dict.subDict("curvatureModel"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::curvatureModel>
Foam::curvatureModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Get subdictionary from the dictionary
    const dictionary coeffDict(dict.subDict("curvatureModel"));

    // Get the name of the desired refinement selection algorithm
    const word curvatureModelTypeName(coeffDict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(curvatureModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "curvatureModel::curvatureModel::New\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Unknown curvatureModel type "
                << curvatureModelTypeName << endl << endl
                << "Valid curvatureModel types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalError);
    }

    return autoPtr<curvatureModel>(cstrIter()(mesh, dict));
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::curvatureModel::~curvatureModel()
{}


// ************************************************************************* //
