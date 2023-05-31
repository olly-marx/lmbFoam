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
    checkProcessorPatchPoints

Description
    Prints out point coordinates on a pair of processors to check sync. Checks
    three things for redundancy:
    1. Points on pointPatch,
    2. Points on polyPatch,
    3. Points on polyPatch by traversing through faces (which will produce a lot
       of duplicates, but here we can see whether the ordering for certain faces
       is messed up).

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "processorPointPatch.H"

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

        // Main part of the name of the files. One for pointPatch and one for
        // polyPatch
        const word baseNamePointPatch
        (
            "pointsOnProcPointPatch" + name(procI) + "-" + name(procJ)
        );

        const word baseNamePolyPatch
        (
            "pointsOnProcPolyPatch" + name(procI) + "-" + name(procJ)
        );

        const word baseNamePolyPatchFaces
        (
            "pointsOnProcPolyPatchFaces" + name(procI) + "-" + name(procJ)
        );

        if (Pstream::myProcNo() == procI || Pstream::myProcNo() == procJ)
        {
            // PART 1. Do the pointPatches first

            // Create a pointMesh and get its boundary
            const pointMesh& pMesh = pointMesh::New(mesh);
            const pointBoundaryMesh& pbMesh = pMesh.boundary();

            // Name of the file
            word procFileName;

            if (Pstream::myProcNo() == procI)
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseNamePointPatch + "-sideOnProc" + name(procI) +".dat");
            }
            else
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseNamePointPatch + "-sideOnProc" + name(procJ) +".dat");
            }

            forAll (pbMesh, patchI)
            {
                const pointPatch& pp = pbMesh[patchI];

                if (isA<processorPointPatch>(pp))
                {
                    const processorPointPatch& ppp =
                        refCast<const processorPointPatch>(pp);

                    // Check whether this is the processor pair that we are
                    // looking for
                    if
                    (
                        (ppp.myProcNo() == procI && ppp.neighbProcNo() == procJ)
                     || (ppp.myProcNo() == procJ && ppp.neighbProcNo() == procI)
                    )
                    {
                        Pout<< "Found processor point patch: " << ppp.name()
                            << nl
                            << "Writing points to file " << procFileName
                            << "..." << endl;

                        // Print out points to the file
                        OFstream of(procFileName);

                        of  << "Points on processor point patch: " << ppp.name()
                            << endl;

                        // Get points on the processor patch
                        const labelList& patchMeshPoints = ppp.meshPoints();

                        // Get points in the mesh
                        const pointField& allPoints = mesh.allPoints();

                        forAll (patchMeshPoints, i)
                        {
                            of  << allPoints[patchMeshPoints[i]] << nl;
                        }

                        of  << endl;
                    }
                } // End if processor point patch
            } // End for all point patches

            // PART 2 and 3: Do the polyPatches next (2 = points directly, 3 =
            // points indirectly through faces)

            // Get polyBoundaryMesh
            const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

            word procFileNameFaces;
            if (Pstream::myProcNo() == procI)
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseNamePolyPatch + "-sideOnProc" + name(procI) +".dat");
                procFileNameFaces = runTime.path()/".."/ // Go out of proc dir
                    (
                        baseNamePolyPatchFaces
                      + "-sideOnProc"
                      + name(procI)
                      + ".dat"
                    );
            }
            else
            {
                procFileName = runTime.path()/".."/ // Go out of proc dir
                    (baseNamePolyPatch + "-sideOnProc" + name(procJ) +".dat");
                procFileNameFaces = runTime.path()/".."/ // Go out of proc dir
                    (
                        baseNamePolyPatchFaces
                      + "-sideOnProc"
                      + name(procJ)
                      + ".dat"
                    );
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
                        Pout<< "Found processor poly patch: " << ppp.name()
                            << nl
                            << "Writing points to file " << procFileName
                            << "..." << endl;

                        // First do the points directly

                        // File to write points to
                        OFstream of(procFileName);

                        of  << "Points on processor poly patch: " << ppp.name()
                            << endl;

                        // Get points on the processor poly patch and write them
                        const pointField& ppPoints = ppp.localPoints();

                        of  << ppPoints << endl;


                        // Now do the points indirectly through faces

                        // File to write points from faces to
                        OFstream off(procFileNameFaces);

                        off << "Face-based points on processor patch: "
                            << ppp.name() << endl;

                        // Grab the points from the mesh
                        const pointField& meshPoints = mesh.allPoints();

                        // Must reverse the face if the processor is slave
                        if (ppp.master())
                        {
                            // Loop through all the faces in the patch
                            forAll (ppp, faceI)
                            {
                                // Get the face
                                const face& f = ppp[faceI];

                                // Loop through the face and print out the
                                // points
                                forAll (f, i)
                                {
                                    off << meshPoints[f[i]] << nl;
                                }

                                off << endl;
                            }
                        }
                        else
                        {
                            // Loop through all the faces in the patch
                            forAll (ppp, faceI)
                            {
                                // Get the reversed face
                                const face rf = ppp[faceI].reverseFace();

                                // Loop through the face and print out the
                                // points
                                forAll (rf, i)
                                {
                                    off << meshPoints[rf[i]] << nl;
                                }

                                off << endl;
                            }
                        }

                        off << endl;
                    }
                } // End if processor point patch
            } // End for all point patches

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
