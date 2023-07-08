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
    Foam::interfaceProximityRefinement

Description
    Selection of refinement cells based on a distance to the interface. Cells
    with cell centres close to the interface get selected for refinement, while
    points away from the interface get selected for unrefinement.

    Relative distance option halves the distance for each current refinement
    level of the cell. Switched on by default.

SourceFiles
    interfaceProximityRefinement.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceProximityRefinement_H
#define interfaceProximityRefinement_H

#include "refinementSelection.H"
#include "interfaceFacesWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class interfaceProximityRefinement Declaration
\*---------------------------------------------------------------------------*/

class interfaceProximityRefinement
:
    public refinementSelection
{
    // Private data

        //- Minimum distance to the interface for refinement in front of the
        //  interface (towards alpha = 0)
        scalar refinementDistanceFront_;

        //- Minimum distance to the interface for refinement back of the
        //  interface (towards alpha = 1)
        scalar refinementDistanceBack_;
        
        //- Minimum distance to the interface for unrefinement in front of the
        //  inteface (towards alpha = 0)
        scalar unrefinementDistanceFront_;

        //- Minimum distance to the interface for unrefinement back of the
        //  inteface (towards alpha = 1)
        scalar unrefinementDistanceBack_;

        //- Relative distance flag: whether to decrease the distance_ for each
        //  cell based on current cell level. Switched on by default
        Switch relativeDistance_;

        //- Interface faces wave object wrapped into autoPtr. Calculated when
        //  calling refinementCellCandidates and deleted when calling
        //  unrefinementPointCandidates (in order to prevent recalculation)
        mutable autoPtr<interfaceFacesWave> interfaceWaveInfo_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceProximityRefinement(const interfaceProximityRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceProximityRefinement&);


public:

    //- Runtime type information
    TypeName("interfaceProximityRefinement");


    // Constructors

        //- Construct from components
        interfaceProximityRefinement(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~interfaceProximityRefinement();


    // Member Functions

        // Selection of refinement/unrefinement candidates

            //- Return transferable list of cells to refine
            virtual Xfer<labelList> refinementCellCandidates() const;

            //- Return transferable list of split points to unrefine
            virtual Xfer<labelList> unrefinementPointCandidates() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
