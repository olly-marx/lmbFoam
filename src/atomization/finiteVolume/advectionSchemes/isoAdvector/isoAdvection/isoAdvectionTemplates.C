/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the IsoAdvector source code library, which is an
    unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "isoAdvection.H"

// ************************************************************************* //

template<typename Type>
Type Foam::isoAdvection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        return f.internalField()[faceI];
    }
    else
    {
        // Boundary face. Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return 0;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "Type isoAdvection::faceValue\n"
                "(\n"
                "    const GeometricField<Type, fvsPatchField, surfaceMesh>&\n"
                "    const label\n"
                ")\n"
            )   << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        const label patchFaceI =
            mesh_.boundaryMesh()[patchI].whichFace(faceI);

        return f.boundaryField()[patchI][patchFaceI];
    }
}


template<typename Type>
void Foam::isoAdvection::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI,
    const Type& value
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        f.internalField()[faceI] = value;
    }
    else
    {
        // Boundary face. Find out which face of which patch
        const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

        // Handle empty patches
        if (mesh_.boundary()[patchI].size() == 0)
        {
            return;
        }

        if (patchI < 0 || patchI >= mesh_.boundaryMesh().size())
        {
            FatalErrorIn
            (
                "void isoAdvection::faceValue\n"
                "(\n"
                "    const GeometricField<Type, fvsPatchField, surfaceMesh>&\n"
                "    const label\n"
                "    const Type\n"
                ")\n"
            )   << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        const label patchFaceI =
            mesh_.boundaryMesh()[patchI].whichFace(faceI);

        f.boundaryField()[patchI][patchFaceI] = value;
    }
}


// ************************************************************************* //
