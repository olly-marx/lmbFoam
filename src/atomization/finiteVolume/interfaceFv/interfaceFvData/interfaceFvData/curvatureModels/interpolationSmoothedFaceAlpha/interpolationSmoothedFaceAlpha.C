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
    Alberto Ceschin, KAUST
\*---------------------------------------------------------------------------*/

#include "interpolationSmoothedFaceAlpha.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interpolationSmoothedFaceAlpha, 0);
addToRunTimeSelectionTable
(
    curvatureModel,
    interpolationSmoothedFaceAlpha,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolationSmoothedFaceAlpha::interpolationSmoothedFaceAlpha
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    curvatureModel(mesh, dict),
    nSmoothIter_
    (
        readLabel(coeffDict().lookup("nSmoothingIters"))
    )
{
    // Sanity check
    if (nSmoothIter_ < 0)
    {
        WarningIn
        (
            "interpolationSmoothedFaceAlpha::interpolationSmoothedFaceAlpha"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict"
            "\n)"
        ) << "Negative number of smoothing iterations specified. "
          << "Disabling smoothing..."
          << endl;

        nSmoothIter_ = 0;
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::interpolationSmoothedFaceAlpha::~interpolationSmoothedFaceAlpha()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interpolationSmoothedFaceAlpha::kappa
(
    const volScalarField& alpha
) const
{
    // Create a copy of the alpha field to prepare for smoothing
    volScalarField alphaSmooth("alphaSmooth", alpha);
    volScalarField alpha_("alphaSmooth",alpha);

    Info<< "Smoothing alpha field " << nSmoothIter_
        << " times for curvature calculation."
        << endl;

    for (label i = 0; i < nSmoothIter_; ++i)
    {
        alphaSmooth = fvc::average(fvc::interpolate(alphaSmooth));
    }
    alphaSmooth.correctBoundaryConditions();

    const volVectorField gradAlpha(fvc::grad(alphaSmooth));

    // Interpolate normals to the faces and scale them
    surfaceVectorField nInterface(fvc::interpolate(gradAlpha));
    nInterface/=
        mag(nInterface) + dimensionedScalar("SMALL", dimless/dimLength, SMALL);

    // Calculate and return the curvature
    tmp<volScalarField> tkappaInterfaceSmooth
    (
        new volScalarField
        (
            "kappaInterface",
            -fvc::div(nInterface & mesh().Sf())
        )
    );

    return tkappaInterfaceSmooth;
}


// ************************************************************************* //
