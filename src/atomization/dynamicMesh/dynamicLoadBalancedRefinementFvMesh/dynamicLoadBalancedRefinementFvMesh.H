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
    Foam::dynamicLoadBalancedRefinementFvMesh

Description
    Dynamic fvMesh with dynamic parallel load balancing and refinement of
    arbitrary polyhedral cells. First, the adaptive refinement is performed,
    followed by optional load balancing, which is done only if the refinement is
    performed. If the refinement has been performed and the ratio of minimum
    number of cells on any processor to maximum number of cells on some other
    processor is smaller than the imbalanceThreshold, the load balancing is
    triggered.

    Derived from dynamicPolyRefinementFvMesh

SourceFiles
    dynamicLoadBalancedRefinementFvMesh.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicLoadBalancedRefinementFvMesh_H
#define dynamicLoadBalancedRefinementFvMesh_H

#include "dynamicPolyRefinementFvMesh.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class dynamicLoadBalancedRefinementFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicLoadBalancedRefinementFvMesh
:
    public dynamicPolyRefinementFvMesh
{

    // Private data

        //- Whether to pefrom parallel load balancing
        Switch loadBalance_;

        //- Imbalance threshold: Minimum allowable ratio between minimum number
        //  of cells on one processor and maximum number of cells cells on the
        //  other processor. If the actual ratio goes below this value, load
        //  balancing is triggered.
        scalar imbalanceThreshold_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicLoadBalancedRefinementFvMesh
        (
            const dynamicLoadBalancedRefinementFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const dynamicLoadBalancedRefinementFvMesh&
        );


public:

    //- Runtime type information
    TypeName("dynamicLoadBalancedRefinementFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicLoadBalancedRefinementFvMesh
        (
            const IOobject& io
        );


    //- Destructor
    virtual ~dynamicLoadBalancedRefinementFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
