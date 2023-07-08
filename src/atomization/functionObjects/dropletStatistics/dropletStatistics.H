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
    dropletStatistics

Description
    Class for gathering droplet statistics.

    Finds disconnected regions of alpha field by:
    1. Marking all cells with alphaThreshold > 1 with increasing index and
       collecting them into a alphaCells dynamic list,
    2. Loop through alphaCells and for each alpha cell, set the alpha index
       of neighbouring cells to the current cell alpha index,
    3. Collect all alpha cells into clusters: all cells with the same alpha
       index (connected cells forming droplets) belong to the same cluster,
    4. For all clusters, calculate the position, velocity and the volume,
    5. Combine cluster data across processor boundaries. The connectivity of all
       clusters on all processors is calculated and the boundary clusters are
       shipped to master processor. Master processor does the agglomeration of
       all clusters spanning across multiple because we might have one cluster
       spanning on more than two processors.

    Wrapped as a function object in dropletStatisticsFunctionObject.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

SourceFiles
    dropletStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef dropletStatistics_H
#define dropletStatistics_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dropletStatistics Declaration
\*---------------------------------------------------------------------------*/

class dropletStatistics
{
    // Helper private typedef
    typedef HashSet<labelPair, Hash<labelPair> > labelPairHashSet;


    // Private data

        //- Name of the function object
        word name_;

        //- Reference to main object registry
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of the alpha field, "alpha1" by default
        word alphaName_;

        //- Threshold alpha value
        scalar alphaThreshold_;

        //- Name of the velocity field, "U" by default
        word UName_;

        //- Minimum droplet diameter to track
        scalar minDiameter_;

        //- Maximum droplet diameter to track
        scalar maxDiameter_;


        // Data members for writing

            //- Directory name to write the data into, "dropletStatistics" by
            //  default
            fileName dirName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dropletStatistics
        (
            const dropletStatistics&
        );

        //- Disallow default bitwise assignment
        void operator=(const dropletStatistics&);


public:

    //- Runtime type information
    TypeName("dropletStatistics");


    // Constructors

        //- Construct from components
        dropletStatistics
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~dropletStatistics();


    // Member Functions

        //- Return name
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, called at each ++ or += of the time-loop
        virtual void execute()
        {
            // Does nothing
        }

        //- Execute at the final time-loop
        virtual void end()
        {
            // Does nothing
        }

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet()
        {
            // Does nothing
        }

        //- Write
        virtual void write();

        //- Update for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {
            // Does nothing
        }

        //- Update for point motion
        virtual void movePoints(const pointField&)
        {
            // Does nothing
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
