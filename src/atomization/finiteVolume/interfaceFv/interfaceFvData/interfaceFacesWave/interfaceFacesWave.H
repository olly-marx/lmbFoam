/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::interfaceFacesWave

Description
    Takes a set of faces defined as interface faces to start MeshWave
    from. After construction, it holds distance at cells.
    The distance can be used to blend free surface sensitive turbulence models
    and redistance the Level Set field.

    Note: in case of multi-region mesh where we don't visit all the cells, we
    will simply set the distance to -1.

SourceFiles
    interfaceFacesWave.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceFacesWave_H
#define interfaceFacesWave_H

#include "cellDistFuncs.H"
#include "FieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wallPoint;
template<class Type> class MeshWave;

/*---------------------------------------------------------------------------*\
                           Class interfaceFacesWave Declaration
\*---------------------------------------------------------------------------*/

class interfaceFacesWave
:
    public cellDistFuncs
{
    // Private Data

        //- Distance at cell centres
        scalarField distance_;


public:

    // Constructors

        //- Construct from mesh. Calculates the distance for all cells
        interfaceFacesWave(const polyMesh& mesh);


    //- Destructor
    virtual ~interfaceFacesWave();


    // Member Functions

        //- Return distance as cell centred field
        const scalarField& distance() const
        {
            return distance_;
        }

        //- Non const access for transferring the contents for efficiency
        scalarField& distance()
        {
            return distance_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
