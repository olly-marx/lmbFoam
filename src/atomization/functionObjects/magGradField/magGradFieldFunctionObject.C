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

#include "magGradFieldFunctionObject.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(magGradFieldFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        magGradFieldFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::magGradFieldFunctionObject::magGradFieldFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(dict.lookupOrDefault<word>("region", polyMesh::defaultRegion)),
    fieldName_(dict.lookup("fieldName")),
    magGradField_
    (
        IOobject
        (
            "magGrad(" + fieldName_ + ")",
            time_.timeName(),
            time_.lookupObject<fvMesh>(regionName_),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        time_.lookupObject<fvMesh>(regionName_),
        dimensionedScalar("zero", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::magGradFieldFunctionObject::~magGradFieldFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::magGradFieldFunctionObject::start()
{
    return true;
}


bool Foam::magGradFieldFunctionObject::execute(const bool)
{
    if (debug)
    {
        Info<< "Updating magGradFieldFunctionObject for field: "
            << fieldName_ << endl;
    }

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    // Search for the field
    if (mesh.foundObject<volScalarField>(fieldName_))
    {
        // Field found, get it
        const volScalarField& field =
            mesh.lookupObject<volScalarField>(fieldName_);

        // Calculate the gradient
        magGradField_ = mag(fvc::grad(field));
        magGradField_.correctBoundaryConditions();

    }
    else
    {
        // Did not find the field, skip
        InfoIn("bool magGradFieldFunctionObject::execute(const bool)")
            << "Field: " << fieldName_ << " not found."
            << nl
            << "Skipping evaluation of gradient..."
            << endl;
    }

    return false;
}


bool Foam::magGradFieldFunctionObject::read(const dictionary& dict)
{
    fieldName_ = word(dict.lookup("fieldName"));

    return false;
}


// ************************************************************************* //
