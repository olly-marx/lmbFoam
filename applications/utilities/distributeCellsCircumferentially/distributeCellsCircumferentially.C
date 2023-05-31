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

Application
    distributeCellsCircumferentially

Description
    A utility that creates a distribution file "circumferentialDistribution" and
    writes there a labelIOList containing information on which cell goes to
    which processor. Cells are collected in circumferential direction in order
    to avoid extreme load imbalance when running AMR cases.

    For example, if one specifies 180 wedges (corresponding to 180 processors),
    cells between 0 and 2 degrees will be distributed to processor 0, 2 to 4
    degrees to processor 1, etc...

Authors
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "labelIOList.H"
#include "cylindricalCS.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Read the dictionary
    IOdictionary distributionDict
    (
        IOobject
        (
            "circumferentialDistributionDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Get axial direction
    vector axialDir = distributionDict.lookup("axialDirection");

    // Sanity check and normalization
    if (mag(axialDir) < SMALL)
    {
        FatalIOErrorIn("distributeCellsCircumferentially", distributionDict)
            << "Zero axialDirection specified. This is not allowed."
            << abort(FatalError);
    }
    else
    {
        axialDir /= mag(axialDir);
    }

    // Get radial direction
    vector radialDir = distributionDict.lookup("radialDirection");

    // Sanity check and normalization
    if (mag(radialDir) < SMALL)
    {
        FatalIOErrorIn("distributeCellsCircumferentially", distributionDict)
            << "Zero radialDirection specified. This is not allowed."
            << abort(FatalError);
    }
    else
    {
        radialDir /= mag(radialDir);
    }

    // Get origin of the coordinate system
    const point origin = distributionDict.lookup("origin");

    // Create the cylindrical coordinate system
    cylindricalCS ccs
    (
        "ccs",
        origin,
        axialDir,
        radialDir
    );

    // Transform cell centres into cylindrical coordinate system
    const vectorField cylindricalCC = ccs.localPosition(mesh.cellCentres());


    // Get number of processors (number of wedges)
    const label nProcs = readLabel(distributionDict.lookup("nProcessors"));

    if (nProcs < 2)
    {
        FatalIOErrorIn("distributeCellsCircumferentially", distributionDict)
            << "Specified less than 2 processors. This is not allowed."
            << abort(FatalError);
    }

    // Calculate delta theta (circumferential sweep) based on number of
    // processors
    const scalar deltaTheta = 360.0/nProcs;

    // Create volScalarField for final decomposition and for visualization
    // purposes
    volScalarField cellCircumferentialDecomposition
    (
        IOobject
        (
            "cellCircumferentialDecomposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("-1", dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );
    scalarField& ccd = cellCircumferentialDecomposition.internalField();

    // Loop through all cell centres in local coordinate system
    forAll (cylindricalCC, cellI)
    {
        // Get theta of this cell. Note: theta defined from -180 to 180, offset
        // it by 180
        const scalar curTheta = cylindricalCC[cellI].y() + 180;

        // Loop through all processors (wedges) to find out where this cell
        // belongs to
        for (label i = 1; i < nProcs + 1; ++i)
        {
            if ((i*deltaTheta >= curTheta) && (curTheta < (i + 1)*deltaTheta))
            {
                // This cell belongs to processor i - 1, mark it and break out
                ccd[cellI] = i - 1;
                break;
            }
        }
    }

    // Sanity checks
    if (max(ccd) > nProcs -1)
    {
        FatalErrorIn("distributeCellsCircumferentially")
            << "Detected maximum processor number: " << max(ccd)
            << " larger than the specified number of processors: " << nProcs
            << abort(FatalError);
    }

    if (min(ccd) == -1)
    {
        FatalErrorIn("distributeCellsCircumferentially")
            << "Found cell that hasn't been distributed to any processor."
            << nl
            << "Check coordinate system definition."
            << abort(FatalError);
    }

    // Convert volScalarField to labelIOList (input for decomposePar with manual
    // decomposition)
    labelIOList circumferentialDistribution
    (
        IOobject
        (
            "circumferentialDistribution",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        ccd.size()
    );

    forAll (ccd, cellI)
    {
        circumferentialDistribution[cellI] = label(ccd[cellI]);
    }

    // Write both the list and the field for visualization. A bit redundant but
    // useful for easy visualization of decomposition for checking
    cellCircumferentialDecomposition.write();
    circumferentialDistribution.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
