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
    Foam::interfaceFvData

Description
    An abstract base class that holds data for implementation of grad(pd) and
    rho jump conditions in the pressure equation, snGrad and grad schemes. The
    class does the following:
    - updates fields needed for blending: wetOwners, wetNeighbours and wetCells,
    - updates interface faces (internal faces dividing submerged and dry cells),
    - updates inverse weighted density and hydrostatic jump fields for interface
      faces.
    - holds old data (previous time step/iteration) that is used by fvc
      (explicit) functions: interfaceSnGrad and grad. This is needed for
      numerical consistency in segregated solvers.
    User of the class has to call member functions that update interface faces
    and weighing factors, update old data and reset old level set field with new
    one for multiple PISO or non orthogonal correctors.
    Level Set and Volume of Fluid variants will be derived from this class.

SourceFiles
    interfaceFvData.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef interfaceFvData_H
#define interfaceFvData_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "MeshObject.H"
#include "mapPolyMesh.H"
#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class interfaceFvData Declaration
\*---------------------------------------------------------------------------*/

class interfaceFvData
:
    public MeshObject<fvMesh, interfaceFvData>
{
public:

    // Public data types

        typedef DynamicList<label> DynamicLabelList;
        typedef PtrList<DynamicLabelList> DynamicLabelListList;


private:

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Model for curvature calculation based on VOF field
        autoPtr<curvatureModel> curvature_;


        // Jump condition members

            //- Density of heavier fluid
            dimensionedScalar rhoPlus_;

            //- Density of lighter fluid
            dimensionedScalar rhoMinus_;

            //- Surface tension between the two fluids
            dimensionedScalar sigma_;

            //- Gravitational acceleration
            vector g_;


            // Current interface data

                //- Surface field denoting if the owner of the face is currently
                //  wet, (1 = wet, 0 = dry)
                surfaceScalarField wetOwners_;

                //- Surface field denoting if the neighbour of the face is
                //  currently wet, (1 = wet, 0 = dry)
                surfaceScalarField wetNeighbours_;

                //- Scalar field for current wet cells (1 = wet, 0 = dry)
                mutable scalarField wetCells_;

                //- DynamicList of labels indicating current interface faces
                mutable DynamicLabelList interfaceFaces_;

                //- DynamicListList of labels indicating current interface faces
                //  at the boundary
                mutable DynamicLabelListList boundaryInterfaceFaces_;

                //- Surface field of current weighted inverse densities
                //  (betaOverbar) for interface faces
                surfaceScalarField betaOverbar_;

                //- Surface field of current hydrostatic jumpss for interface
                //  faces
                surfaceScalarField hydrostaticJump_;

                //- Position of the interface for this face. Needed to calculate
                //  the distance from the interface using interfaceFacesWave
                surfaceVectorField cGamma_;


            // Old interface data

                //- Surface field denoting if the owner of the face was wet in
                //  the previous time step/iteration, (1 = wet, 0 = dry)
                surfaceScalarField wetOwnersOld_;

                //- Volume field for old wet cells (1 = wet, 0 = dry)
                mutable scalarField wetCellsOld_;

                //- DynamicList of labels indicating old interface faces
                mutable DynamicLabelList interfaceFacesOld_;

                //- DynamicListList of labels indicating old interface faces
                //  at the boundary
                mutable DynamicLabelListList boundaryInterfaceFacesOld_;

                //- Surface field of old weighted inverse densities
                //  (betaOverbar) for interface faces
                surfaceScalarField betaOverbarOld_;

                //- Surface field of old hydrostatic jumps for interface faces
                surfaceScalarField hydrostaticJumpOld_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            interfaceFvData(const interfaceFvData&);

            //- Disallow default bitwise assignment
            void operator=(const interfaceFvData&);


protected:

    // Protected Member Functions

        //- Non-const access to the curvature model
        inline curvatureModel& curvature();


        // Access to current interface finite volume data

            //- Return access to current wet owners
            inline surfaceScalarField& wetOwners();

            //- Return access to current wet neighbours
            inline surfaceScalarField& wetNeighbours();

            //- Return access to current wet cells
            inline scalarField& wetCells();

            //- Return access to current interface faces
            inline DynamicLabelList& interfaceFaces();

            //- Return access to current boundary interface faces
            inline DynamicLabelListList& boundaryInterfaceFaces();

            //- Return access to current weighted inverse density field
            inline surfaceScalarField& betaOverbar();

            //- Return access to current hydrostaticJump field
            inline surfaceScalarField& hydrostaticJump();

            //- Return access to current cGamma field
            inline surfaceVectorField& cGamma();


        // Access to old interface finite volume data

            //- Return access to old wet owners
            inline surfaceScalarField& wetOwnersOld();

            //- Return access to current wet cells
            inline scalarField& wetCellsOld();

            //- Return access to old interface faces
            inline DynamicLabelList& interfaceFacesOld();

            //- Return access to old boundary interface faces
            inline DynamicLabelListList& boundaryInterfaceFacesOld();

            //- Return access to old weighted inverse density field
            inline surfaceScalarField& betaOverbarOld();

            //- Return access to old hydrostaticJump field
            inline surfaceScalarField& hydrostaticJumpOld();


public:

    // Declare name of the class and its debug switch
    TypeName("interfaceFvData")


    //- Constructor (from mesh)
    interfaceFvData(const fvMesh& mesh);


    //- Destructor
    virtual ~interfaceFvData()
    {}


    // Member Functions

        // Checking helper functions

            //- Check whether interface corrected schemes are specified in
            //  fvSchemes
            void checkSchemes
            (
                const word& laplacianSchemeName,
                const word& snGradSchemeName,
                const word& gradSchemeName,
                const bool ghostFluidMethod
            ) const;


        // Access

            //- Return the reference to mesh
            inline const fvMesh& mesh() const;
            
            //- Const access to the curvature model
            inline const curvatureModel& curvature() const;

            //- Return density of heavier fluid
            inline const dimensionedScalar& rhoPlus() const;

            //- Return density of lighter fluid
            inline const dimensionedScalar& rhoMinus() const;

            //- Return inverse density of heavier fluid
            inline scalar betaPlus() const;

            //- Return inverse density of lighter fluid
            inline scalar betaMinus() const;

            //- Return surface tension coefficient between the fluids
            inline const dimensionedScalar& sigma() const;

            //- Return gravitational acceleration vector
            inline vector g() const;


            // Access to current interface finite volume data

                //- Return const access to current wet owners
                inline const surfaceScalarField& wetOwners() const;

                //- Return const access to current wet neighbours
                inline const surfaceScalarField& wetNeighbours() const;

                //- Return const access to current wet cells
                inline const scalarField& wetCells() const;

                //- Return const access to current interface faces
                inline const DynamicLabelList& interfaceFaces() const;

                //- Return const access to current boundary interface faces
                inline const DynamicLabelListList&
                boundaryInterfaceFaces() const;

                //- Return const access to current weighted inverse density
                //  field
                inline const surfaceScalarField& betaOverbar() const;

                //- Return const access to current hydrostaticJump field
                inline const surfaceScalarField& hydrostaticJump() const;

                //- Return const access to current cGamma field
                inline const surfaceVectorField& cGamma() const;


            // Access to old interface finite volume data

                //- Return const access to old wet owners
                inline const surfaceScalarField& wetOwnersOld() const;

                //- Return const access to current wet cells
                inline const scalarField& wetCellsOld() const;

                //- Return const access to old interface faces
                inline const DynamicLabelList& interfaceFacesOld() const;

                //- Return const access to old boundary interface faces
                inline const DynamicLabelListList&
                boundaryInterfaceFacesOld() const;

                //- Return const access to old weighted inverse density field
                inline const surfaceScalarField& betaOverbarOld() const;

                //- Return const access to old hydrostaticJump field
                inline const surfaceScalarField& hydrostaticJumpOld() const;


        // Jump condition corrections

            //- Update all blending fields
            virtual void update() = 0;

            //- Update old interface data to correspond to current state
            virtual void updateOld();


        // Mesh change corrections

            //- Mesh motion update. Does nothing
            virtual bool movePoints() const
            {
                return false;
            }

            //- Topo change update. Resets the size of necessary fields
            //  (e.g. boundaryInterfaceFaces need to be updated in DLB runs)
            virtual bool updateMesh(const mapPolyMesh& map) const;


        //- Write data. Does nothing
        virtual bool writeData(Ostream&) const
        {
            return false;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceFvDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
