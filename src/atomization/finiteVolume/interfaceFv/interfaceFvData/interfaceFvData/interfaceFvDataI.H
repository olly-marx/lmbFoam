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

\*---------------------------------------------------------------------------*/

#include "interfaceFvData.H"

// * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //

inline Foam::curvatureModel& Foam::interfaceFvData::curvature()
{
    return curvature_();
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::wetOwners()
{
    return wetOwners_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::wetNeighbours()
{
    return wetNeighbours_;
}


inline Foam::scalarField& Foam::interfaceFvData::wetCells()
{
    return wetCells_;
}


inline Foam::interfaceFvData::DynamicLabelList&
Foam::interfaceFvData::interfaceFaces()
{
    return interfaceFaces_;
}


inline Foam::interfaceFvData::DynamicLabelListList&
Foam::interfaceFvData::boundaryInterfaceFaces()
{
    return boundaryInterfaceFaces_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::betaOverbar()
{
    return betaOverbar_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::hydrostaticJump()
{
    return hydrostaticJump_;
}


inline Foam::surfaceVectorField& Foam::interfaceFvData::cGamma()
{
    return cGamma_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::wetOwnersOld()
{
    return wetOwnersOld_;
}


inline Foam::scalarField& Foam::interfaceFvData::wetCellsOld()
{
    return wetCellsOld_;
}


inline Foam::interfaceFvData::DynamicLabelList&
Foam::interfaceFvData::interfaceFacesOld()
{
    return interfaceFacesOld_;
}


inline Foam::interfaceFvData::DynamicLabelListList&
Foam::interfaceFvData::boundaryInterfaceFacesOld()
{
    return boundaryInterfaceFacesOld_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::betaOverbarOld()
{
    return betaOverbarOld_;
}


inline Foam::surfaceScalarField& Foam::interfaceFvData::hydrostaticJumpOld()
{
    return hydrostaticJumpOld_;
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::interfaceFvData::mesh() const
{
    return mesh_;
}


inline const Foam::curvatureModel& Foam::interfaceFvData::curvature() const
{
    return curvature_();
}


inline const Foam::dimensionedScalar& Foam::interfaceFvData::rhoPlus() const
{
    return rhoPlus_;
}


inline const Foam::dimensionedScalar& Foam::interfaceFvData::rhoMinus() const
{
    return rhoMinus_;
}


inline Foam::scalar Foam::interfaceFvData::betaPlus() const
{
    return 1/rhoPlus_.value();
}


inline Foam::scalar Foam::interfaceFvData::betaMinus() const
{
    return 1/rhoMinus_.value();
}


inline const Foam::dimensionedScalar& Foam::interfaceFvData::sigma() const
{
    return sigma_;
}


inline Foam::vector Foam::interfaceFvData::g() const
{
    return g_;
}


inline const Foam::surfaceScalarField& Foam::interfaceFvData::wetOwners() const
{
    return wetOwners_;
}


inline const Foam::surfaceScalarField& Foam::interfaceFvData::wetNeighbours()
const
{
    return wetNeighbours_;
}


inline const Foam::scalarField& Foam::interfaceFvData::wetCells() const
{
    return wetCells_;
}


inline const Foam::interfaceFvData::DynamicLabelList&
Foam::interfaceFvData::interfaceFaces() const
{
    return interfaceFaces_;
}


inline const Foam::interfaceFvData::DynamicLabelListList&
Foam::interfaceFvData::boundaryInterfaceFaces() const
{
    return boundaryInterfaceFaces_;
}


inline const Foam::surfaceScalarField&
Foam::interfaceFvData::betaOverbar() const
{
    return betaOverbar_;
}


inline const Foam::surfaceScalarField&
Foam::interfaceFvData::hydrostaticJump() const
{
    return hydrostaticJump_;
}


inline const Foam::surfaceVectorField&
Foam::interfaceFvData::cGamma() const
{
    return cGamma_;
}


inline const Foam::surfaceScalarField&
Foam::interfaceFvData::wetOwnersOld() const
{
    return wetOwnersOld_;
}


inline const Foam::scalarField& Foam::interfaceFvData::wetCellsOld() const
{
    return wetCellsOld_;
}


inline const Foam::interfaceFvData::DynamicLabelList&
Foam::interfaceFvData::interfaceFacesOld() const
{
    return interfaceFacesOld_;
}


inline const Foam::interfaceFvData::DynamicLabelListList&
Foam::interfaceFvData::boundaryInterfaceFacesOld() const
{
    return boundaryInterfaceFacesOld_;
}


inline const Foam::surfaceScalarField&
Foam::interfaceFvData::betaOverbarOld() const
{
    return betaOverbarOld_;
}


inline const Foam::surfaceScalarField&
Foam::interfaceFvData::hydrostaticJumpOld() const
{
    return hydrostaticJumpOld_;
}


// ************************************************************************* //
