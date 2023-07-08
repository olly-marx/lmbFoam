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
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "interpolationSmoothedLevelSet.H"
#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "interfaceFacesWave.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "pointVolInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interpolationSmoothedLevelSet, 0);
addToRunTimeSelectionTable
(
    curvatureModel,
    interpolationSmoothedLevelSet,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolationSmoothedLevelSet::interpolationSmoothedLevelSet
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
            "interpolationSmoothedLevelSet::interpolationSmoothedLevelSet"
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

Foam::interpolationSmoothedLevelSet::~interpolationSmoothedLevelSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interpolationSmoothedLevelSet::kappa
(
    const volScalarField& alpha
) const
{
    // Create interface faces object and get the distance to the interface
    const interfaceFacesWave interfaceWaveInfo(mesh());
    const scalarField& interfaceDistance = interfaceWaveInfo.distance();

    // Calculate level set field from distance field and alpha. Note: we use
    // zeroGradient BCs all around for the first iteration. We'll then calculate
    // the gradient of levelSet with which we'll correct the boundary values in
    // a new field with calculated BCs
    volScalarField levelSetZG
    (
        IOobject
        (
            "levelSetZG",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    levelSetZG.internalField() =
        pos(alpha - 0.5)*interfaceDistance - neg(alpha - 0.5)*interfaceDistance;
    levelSetZG.correctBoundaryConditions();

    // Calculate the gradient of level set field
    volVectorField gradls("gradlsZG", fvc::grad(levelSetZG));
    gradls.write();

    // Create a new field with different BCs
    volScalarField levelSet
    (
        IOobject
        (
            "levelSet",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        levelSetZG,
        calculatedFvPatchScalarField::typeName
    );

    // Loop through all the non-coupled boundaries
    auto& blevelSet = levelSet.boundaryField();
    const auto& bgradls = gradls.boundaryField();
    forAll (blevelSet, patchI)
    {
        // Get level set field at the patch
        fvPatchScalarField& pls = blevelSet[patchI];

        // Check whether the patch is coupled
        if (!pls.coupled())
        {
            // Get gradient and delta vectors for extrapolation
            const vectorField pgradls = bgradls[patchI].patchInternalField();
            const vectorField pdelta = blevelSet[patchI].patch().delta();

            // Extrapolate from cell centres to patches. Note: patch level set
            // field (pls) already holds cell centred values because it was
            // created with zeroGradient BCs. Hence the += operator
            pls += pdelta & pgradls;
        }
    }

    // Smooth out level set by interpolating onto points and then back to cell
    // centres a number of times
    if (nSmoothIter_ != 0)
    {
        Info<< "Smoothing level set field " << nSmoothIter_
            << " times for curvature calculation."
            << endl;

        // Create volume to point and point to volume interpolation objects
        const volPointInterpolation& vpi = volPointInterpolation::New(mesh());
        const pointMesh& pMesh = pointMesh::New(mesh());
        const pointVolInterpolation pvi(pMesh, mesh());

        for (label i = 0; i < nSmoothIter_; ++i)
        {
            // Interpolate to points and then back to cell centres
            pointScalarField levelSetp(vpi.interpolate(levelSet));
            levelSet = pvi.interpolate(levelSetp);
        }
    }

    // Now that we have the level set field with correct boundaries,
    // re-calculate the gradient
    gradls = fvc::grad(levelSet);

    // Calculate the normals to the faces and scale them
    surfaceVectorField nInterface(fvc::interpolate(gradls));
    nInterface/=
        mag(nInterface) + dimensionedScalar("SMALL", dimless/dimLength, SMALL);

    // Calculate and return the curvature
    tmp<volScalarField> tkappaInterface
    (
        new volScalarField
        (
            "kappaInterface",
            -fvc::div(nInterface & mesh().Sf())
        )
    );

    return tkappaInterface;
}


// ************************************************************************* //
