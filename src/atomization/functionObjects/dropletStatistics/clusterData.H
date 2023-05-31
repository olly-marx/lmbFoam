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
    Foam::clusterData

Description
    A class that holds cluster data: volume of the cluster, centre of
    volume, velocity and effective radius (radius assuming that the cluster
    forms a spherical droplet).

SourceFiles

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef clusterData_H
#define clusterData_H

#include "label.H"
#include "point.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;
class clusterData;

// Forward declaration of friend operators
inline bool operator==(const clusterData&, const clusterData&);
inline bool operator!=(const clusterData&, const clusterData&);

Istream& operator>>(Istream&, clusterData&);
Ostream& operator<<(Ostream&, const clusterData&);


/*---------------------------------------------------------------------------*\
                        Class clusterData Declaration
\*---------------------------------------------------------------------------*/

class clusterData
{
    // Private data

        // Cluster characteristics

            //- Volume of the cluster
            scalar volume_;

            //- Centre of volume
            point centre_;

            //- Velocity of the cluster
            vector velocity_;

            //- Flag whether the cluster is active (for easier agglomeration
            //  across processor patches)
            bool active_;


public:

    // Constructors

        //- Construct null
        inline clusterData()
        :
            volume_(-GREAT),
            centre_(point::min),
            velocity_(vector::zero),
            active_(false)
        {}

        //- Construct from components
        inline clusterData
        (
            const scalar volume,
            const point& centre,
            const vector& velocity,
            const bool active = true
        )
        :
            volume_(volume),
            centre_(centre),
            velocity_(velocity),
            active_(active)
        {}

        //- Construct from IStream
        inline clusterData(Istream& is)
        :
            volume_(readLabel(is)),
            centre_(is),
            velocity_(is),
            active_(readBool(is))
        {}

        //- Copy constructor: default


    //- Destructor: default


    // Member Functions

        // Access functions

            //- Return const access to volume
            const scalar& volume() const
            {
                return volume_;
            }

            //- Return access to volume
            scalar& volume()
            {
                return volume_;
            }

            //- Return const access to centre
            const point& centre() const
            {
                return centre_;
            }

            //- Return access to centre
            point& centre()
            {
                return centre_;
            }

            //- Return const access to velocity
            const vector& velocity() const
            {
                return velocity_;
            }

            //- Return access to velocity
            vector& velocity()
            {
                return velocity_;
            }

            //- Check whether the cluster is active
            bool active() const
            {
                return active_;
            }


        //- Calculate and return effective diameter
        scalar effectiveDiameter() const
        {
            return Foam::pow(6.0*volume_/mathematicalConstant::pi, 1.0/3.0);
        }

        //- Whether the diameter of this cluster is within given bounds
        bool diameterWithinBounds(const scalar& dMin, const scalar& dMax) const
        {
            const scalar d = effectiveDiameter();

            return (d > dMin) && (d < dMax);
        }


        // Edit

            //- Set data
            void setData
            (
                const scalar volume,
                const point& centre,
                const vector& velocity,
                const bool active = true
            )
            {
                volume_ = volume;
                centre_ = centre;
                velocity_ = velocity;
                active_ = active;
            }

            //- Deactivate
            void deactivate()
            {
                active_ = false;
            }


        // Friend operators

            //- Equal comparison operator. Note: needed to allow instantiation
            //  of IOstream operators for lists. No two clusters are alike.
            friend bool operator==(const clusterData& a, const clusterData& b)
            {
                return false;
            }

            //- Not-equal comparison operator. Note: needed to allow
            //  instantiation of IOstream operators for lists. Implemeted in
            //  terms of operator==
            friend bool operator!=(const clusterData& a, const clusterData& b)
            {
                return !(a == b);
            }


        // IOstream operators

            //- Istream operator
            friend Istream& operator>>(Istream& is, clusterData& cd)
            {
                is  >> cd.volume_
                    >> cd.centre_
                    >> cd.velocity_
                    >> cd.active_;

                return is;
            }

            //- Ostream operator
            friend Ostream& operator<<(Ostream& os, const clusterData& cd)
            {
                os  << cd.volume_ << token::SPACE
                    << cd.centre_ << token::SPACE
                    << cd.velocity_ << token::SPACE
                    << cd.active_;

                return os;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
