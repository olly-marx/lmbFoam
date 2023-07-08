/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::magGradFieldFunctionObject

SourceFiles
    magGradFieldFunctionObject.C

Description
    Function object responsible for creating a volScalarField mag(grad(field))
    on the fly.

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef magGradFieldFunctionObject_H
#define magGradFieldFunctionObject_H

#include "functionObject.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                 Class magGradFieldFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class magGradFieldFunctionObject
:
    public functionObject
{
    // Private data

        //- Name of this set of magGradFieldFunctionObject objects
        const word name_;

        //- Object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Name of the field for gradient calculation
        word fieldName_;

        //- Gradient field
        volScalarField magGradField_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        magGradFieldFunctionObject(const magGradFieldFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const magGradFieldFunctionObject&);


public:

    //- Runtime type information
    TypeName("magGradField");


    // Constructors

        //- Construct from compnents
        magGradFieldFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Destructor

        virtual ~magGradFieldFunctionObject();


    // Member Functions

        //- Name
        virtual const word& name() const
        {
            return name_;
        }

        //- Start
        virtual bool start();

        //- Execute
        virtual bool execute(const bool forceWrite);

        //- Read the magGradFieldFunctionObject data
        virtual bool read(const dictionary& dict);

        //- Update for topo mesh changes
        virtual void updateMesh(const mapPolyMesh& mpm)
        {
            // Does nothing

            // Note: magGradField_ is already updates when fvMesh::updateMesh is
            // called since it is in the registry
        }

        //- Update for changes of point motion
        virtual void movePoints(const pointField& points)
        {
            // Does nothing
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
