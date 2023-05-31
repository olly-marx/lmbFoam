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

#include "interpolationSmoothedCSFAlpha.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interpolationSmoothedCSFAlpha, 0);
addToRunTimeSelectionTable
(
    curvatureModel,
    interpolationSmoothedCSFAlpha,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolationSmoothedCSFAlpha::interpolationSmoothedCSFAlpha
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
        WarningInFunction
            << "Negative number of smoothing iterations specified. "
            << "Disabling smoothing..."
            << endl;

        nSmoothIter_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interpolationSmoothedCSFAlpha::kappa
(
    const volScalarField& alpha
) const
{
    volScalarField alphaSmooth("alphaSmooth", alpha);
    volScalarField alphaUnderscore("alphaSmooth", alpha);

    Info<< "Smoothing alpha field " << nSmoothIter_
        << " times for curvature calculation."
        << endl;

    const scalar cSK = 0.5;

    for (label i = 0; i < nSmoothIter_; ++i)
    {
        alphaSmooth =
            cSK*(fvc::average(fvc::interpolate(alphaSmooth)))
          + (scalar(1) - cSK)*alphaSmooth;
    }

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
    volScalarField& kappaInterfaceSmooth = tkappaInterfaceSmooth();

    volScalarField kappaInterfaceUnderscore =
        -fvc::div(nInterface & mesh().Sf());
    
    volScalarField w =
        sqrt(alphaUnderscore*(scalar(1) - alphaUnderscore)
      + dimensionedScalar("SMALL", dimless, SMALL));
    
    volScalarField wSmooth = fvc::average(fvc::interpolate(w));

    for (label i = 0; i < nSmoothIter_; ++i)
    {    
        kappaInterfaceSmooth =
            (fvc::average(fvc::interpolate(kappaInterfaceSmooth*w)))/wSmooth;

        kappaInterfaceSmooth =
            2*sqrt(alphaUnderscore*(1-alphaUnderscore))*
            kappaInterfaceUnderscore
          + (1 - 2*sqrt(alphaUnderscore*(1.-alphaUnderscore)))*
            kappaInterfaceSmooth;
   }

    // Correct coupled patches.  HJ, 14/Oct/2021
    kappaInterfaceSmooth.correctBoundaryConditions();
    
    return tkappaInterfaceSmooth;
}


// ************************************************************************* //
