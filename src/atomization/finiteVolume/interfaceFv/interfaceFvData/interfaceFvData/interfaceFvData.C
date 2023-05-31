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
#include "uniformDimensionedFields.H"
#include "interfaceGaussLaplacian.H"
#include "interfaceSnGradScheme.H"
#include "interfaceGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interfaceFvData, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceFvData::interfaceFvData(const fvMesh& mesh)
:
    MeshObject<fvMesh, interfaceFvData>(mesh),

    mesh_(mesh),
    curvature_(), // Initialized in constructor body

    rhoPlus_(dimensionedScalar("zero", dimDensity, 0)),
    rhoMinus_(dimensionedScalar("zero", dimDensity, 0)),
    sigma_(dimensionedScalar("zero", dimMass/dimTime/dimTime, 0)),
    g_
    (
        mesh.thisDb().lookupObject<uniformDimensionedVectorField>
        (
            "g"
        ).value()
    ),

    wetOwners_
    (
        IOobject
        (
            "wetOwners",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    wetNeighbours_
    (
        IOobject
        (
            "wetNeighbours",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    wetCells_(mesh.nCells()),
    // Guess the size to prevent excessive reallocation
    interfaceFaces_(label(pow(mesh.nFaces(), 2.0/3.0))),
    boundaryInterfaceFaces_(mesh.boundary().size()),
    betaOverbar_
    (
        IOobject
        (
            "betaOverbar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    hydrostaticJump_
    (
        IOobject
        (
            "hydrostaticJump",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    cGamma_
    (
        IOobject
        (
            "cGamma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    ),

    wetOwnersOld_
    (
        IOobject
        (
            "wetOwnersOld",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    wetCellsOld_(mesh.nCells()),
    // Guess the size to prevent excessive reallocation
    interfaceFacesOld_(label(pow(mesh.nFaces(), 2.0/3.0))),
    boundaryInterfaceFacesOld_(mesh.boundary().size()),
    betaOverbarOld_
    (
        IOobject
        (
            "betaOverbarOld",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    ),
    hydrostaticJumpOld_
    (
        IOobject
        (
            "hydrostaticJumpOld",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    )
{
    // Read densities from transport properties
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    rhoPlus_ = dimensionedScalar
    (
        transportProperties.subDict("phase1").lookup("rho")
    );

    rhoMinus_ = dimensionedScalar
    (
        transportProperties.subDict("phase2").lookup("rho")
    );

    sigma_ = dimensionedScalar
    (
        transportProperties.lookup("sigma")
    );

    // Set boundary interface faces. Since we do not expect a lot of interface
    // faces at the boundaries, size of DynamicLabelLists for patches is
    // initialized with 128
    forAll(boundaryInterfaceFaces_, patchI)
    {
        boundaryInterfaceFaces_.set(patchI, new DynamicLabelList(128));
        boundaryInterfaceFacesOld_.set(patchI, new DynamicLabelList(128));
    }

    // Create the curvature model
    curvature_ = curvatureModel::New(mesh, transportProperties);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interfaceFvData::checkSchemes
(
    const word& laplacianSchemeName,
    const word& snGradSchemeName,
    const word& gradSchemeName,
    const bool ghostFluidMethod
) const
{
    if (ghostFluidMethod)
    {
        // Check if the interface corrected scheme is used for the laplacian of
        // pressure
        const tmp<fv::laplacianScheme<scalar, scalar> > pdLaplacianScheme
        (
            fv::laplacianScheme<scalar, scalar>::New
            (
                mesh_,
                mesh_.schemesDict().laplacianScheme(laplacianSchemeName)
            )
        );

        if (!isA<fv::interfaceGaussLaplacian>(pdLaplacianScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "Laplacian scheme " << laplacianSchemeName
              << " is not interface corrected." << nl
              << "Please use interface corrected laplacian scheme "
              << "in order to use the Ghost Fluid Method." << endl
              << exit(FatalError);
        }

        // Check if the interface corrected scheme is used for the snGrad of
        // pressure
        const tmp<fv::snGradScheme<scalar> > pdSnGradScheme
        (
            fv::snGradScheme<scalar>::New
            (
                mesh_,
                mesh_.schemesDict().snGradScheme(snGradSchemeName)
            )
        );

        if (!isA<fv::interfaceSnGradScheme>(pdSnGradScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "snGrad scheme " << snGradSchemeName
              << " is not interface corrected." << nl
              << "Please use interface corrected snGrad scheme in order to use "
              << "the Ghost Fluid Method." << endl
              << exit(FatalError);
        }

        // Check if the interface gradient scheme is used for the grad of
        // pressure
        const tmp<fv::gradScheme<scalar> > pdGradScheme
        (
            fv::gradScheme<scalar>::New
            (
                mesh_,
                mesh_.schemesDict().gradScheme(gradSchemeName)
            )
        );

        if (!isA<fv::interfaceGradScheme>(pdGradScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "grad scheme " << gradSchemeName
              << " is not interface corrected." << nl
              << "Please use interface corrected grad scheme in order to use "
              << "the Ghost Fluid Method." << endl
              << exit(FatalError);
        }
    }
    else
    {
        // Check if the ordinary scheme is used for the laplacian of
        // pressure
        const tmp<fv::laplacianScheme<scalar, scalar> > pdLaplacianScheme
        (
            fv::laplacianScheme<scalar, scalar>::New
            (
                mesh_,
                mesh_.schemesDict().laplacianScheme(laplacianSchemeName)
            )
        );

        if (isA<fv::interfaceGaussLaplacian>(pdLaplacianScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "Laplacian scheme " << laplacianSchemeName
              << " is interface corrected." << nl
              << "Please use ordinary laplacian scheme "
              << "in case you are not using the Ghost Fluid Method." << endl
              << exit(FatalError);
        }

        // Check if the ordinary snGrad scheme is used for the snGrad of
        // pressure
        const tmp<fv::snGradScheme<scalar> > pdSnGradScheme
        (
            fv::snGradScheme<scalar>::New
            (
                mesh_,
                mesh_.schemesDict().snGradScheme(snGradSchemeName)
            )
        );

        if (isA<fv::interfaceSnGradScheme>(pdSnGradScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "snGrad scheme " << snGradSchemeName
              << " is interface corrected." << nl
              << "Please use ordinary snGrad scheme in case you are not using "
              << "the Ghost Fluid Method." << endl
              << exit(FatalError);
        }

        // Check if the ordinary gradient scheme is used for the grad of
        // pressure
        const tmp<fv::gradScheme<scalar> > pdGradScheme
        (
            fv::gradScheme<scalar>::New
            (
                mesh_,
                mesh_.schemesDict().gradScheme(gradSchemeName)
            )
        );

        if (isA<fv::interfaceGradScheme>(pdGradScheme()))
        {
            FatalErrorIn
            (
                "Foam::interfaceFvData::checkSchemes"
                "(\n"
                "    const word& laplacianSchemeName,\n"
                "    const word& snGradSchemeName,\n"
                "    const word& gradSchemeName,\n"
                "    const bool ghostFluidMethod\n"
                ")\n const"
            ) << "grad scheme " << gradSchemeName
              << " is interface corrected." << nl
              << "Please use ordinary grad scheme in case you are not using "
              << "the Ghost Fluid Method." << endl
              << exit(FatalError);
        }
    }
}


void Foam::interfaceFvData::updateOld()
{
    wetCellsOld_ = wetCells_;

    wetOwnersOld_ == wetOwners_;
    betaOverbarOld_ == betaOverbar_;
    hydrostaticJumpOld_ == hydrostaticJump_;

    interfaceFacesOld_ = interfaceFaces_;
    forAll(boundaryInterfaceFacesOld_, patchI)
    {
        boundaryInterfaceFacesOld_[patchI] = boundaryInterfaceFaces_[patchI];
    }
}


bool Foam::interfaceFvData::updateMesh(const mapPolyMesh& map) const
{
    // Get number of faces to reserve. Note: conversion to label
    const label reservedSize = pow(mesh_.nFaces(), 2.0/3.0);

    // Clear the list and reserve enough storage
    interfaceFaces_.clear();
    interfaceFaces_.reserve(reservedSize);

    interfaceFacesOld_.clear();
    interfaceFacesOld_.reserve(reservedSize);

    // Clear all boundary interface faces and set size. Note: needed for
    // parallel load balancing runs where number of patches might change
    boundaryInterfaceFaces_.clear();
    boundaryInterfaceFaces_.setSize(map.mesh().boundaryMesh().size());

    boundaryInterfaceFacesOld_.clear();
    boundaryInterfaceFacesOld_.setSize(map.mesh().boundaryMesh().size());

    forAll (boundaryInterfaceFaces_, patchI)
    {
        boundaryInterfaceFaces_.set(patchI, new DynamicLabelList(128));
        boundaryInterfaceFacesOld_.set(patchI, new DynamicLabelList(129));
    }

    // Update wetCells as well (ordinary scalarField)
    wetCells_.setSize(map.mesh().nCells());
    wetCellsOld_.setSize(map.mesh().nCells());

    return true;
}


// ************************************************************************* //
