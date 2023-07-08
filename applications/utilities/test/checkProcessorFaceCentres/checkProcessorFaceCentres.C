/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Application
    checkProcessorFaceCentres

Description
    Prints out face centres on a pair of processors to check sync.

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert("procI", "label");
    argList::validOptions.insert("procJ", "label");

#   include "setRootCase.H"
#   include "createTime.H"

    if (!Pstream::parRun())
    {
        FatalErrorIn(args.executable())
            << "It only makes sense to run procIndex for a parallel run."
            << exit(FatalError);
    }

#   include "createMesh.H"

    if (args.optionFound("procI") && args.optionFound("procJ"))
    {
        const label procI = readLabel(IStringStream(args.options()["procI"])());
        const label procJ = readLabel(IStringStream(args.options()["procJ"])());

        // Main part of the name of the file
        const word baseName
        (
            "faceCentresOnProcPatch" + name(procI) + "-" + name(procJ)
        );

        if (Pstream::myProcNo() == procI || Pstream::myProcNo() == procJ)
        {
            // Get boundary mesh
            const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

            // Name of the file
            word procFileName;

            if (Pstream::myProcNo() == procI)
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseName + "-sideOnProc" + name(procI) +".dat");
            }
            else
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseName + "-sideOnProc" + name(procJ) +".dat");
            }

            forAll (bMesh, patchI)
            {
                const polyPatch& pp = bMesh[patchI];

                if (isA<processorPolyPatch>(pp))
                {
                    const processorPolyPatch& ppp =
                        refCast<const processorPolyPatch>(pp);

                    // Check whether this is the processor pair that we are
                    // looking for
                    if
                    (
                        (ppp.myProcNo() == procI && ppp.neighbProcNo() == procJ)
                     || (ppp.myProcNo() == procJ && ppp.neighbProcNo() == procI)
                    )
                    {
                        Pout<< "Found processor patch: " << ppp.name()
                            << nl
                            << "Writing points to file " << procFileName
                            << "..." << endl;

                        // Print out points to the file
                        OFstream of(procFileName);

                        of  << "Face centres on processor patch: " << ppp.name()
                            << nl << ppp.faceCentres()
                            << endl;
                    }
                } // End if processor patch
            } // End for all patches
        } // End if these are processors we're looking for
    }
    else
    {
        FatalErrorIn(args.executable())
            << "You must provide processor pair to test: procI and procJ."
            << nl
            << "Example usage: "
            << "mpirun -np 4 checkProcessorPatchPoints "
            << "-procI 1 -procJ 3 -parallel"
            << exit(FatalError);
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
