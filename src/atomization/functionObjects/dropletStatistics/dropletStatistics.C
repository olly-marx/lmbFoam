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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "dropletStatistics.H"
#include "volFields.H"
#include "clusterData.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"
#include "processorFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dropletStatistics, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dropletStatistics::dropletStatistics
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),

    alphaName_(dict.lookupOrDefault<word>("alphaName", "alpha1")),
    alphaThreshold_(readScalar(dict.lookup("alphaThreshold"))),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    minDiameter_(readScalar(dict.lookup("minDiameter"))),
    maxDiameter_(readScalar(dict.lookup("maxDiameter"))),

    dirName_(dict.lookupOrDefault<word>("dirName", "dropletStatistics"))
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "dropletStatistics::dropletStatistics"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

    Info<< "Creating dropletStatistics function object. " << endl;

    // Create directory dropletStatistics directory if not already present
    fileName outputDir;

    if (Pstream::parRun())
    {
        outputDir = obr.time().rootPath()/obr.time().caseName()/".."/dirName_;
    }
    else
    {
        outputDir = obr.time().rootPath()/obr.time().caseName()/dirName_;
    }

    if (!exists(outputDir))
    {
        mkDir(outputDir);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dropletStatistics::~dropletStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dropletStatistics::read(const dictionary& dict)
{
    // Read updates if necessary
    if (active_)
    {
        alphaName_ = dict.lookupOrDefault<word>("alphaName", "alpha1");
        alphaThreshold_ = dict.lookupOrDefault<scalar>("alphaThreshold", 0.1);
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        minDiameter_ = readScalar(dict.lookup("minDiameter"));
        maxDiameter_ = readScalar(dict.lookup("maxDiameter"));
        dirName_ = dict.lookupOrDefault<word>("dirName", "dropletStatistics");
    }
}


void Foam::dropletStatistics::write()
{
    if (!active_)
    {
        return;
    }

    // Get fvMesh
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if
    (
        mesh.foundObject<volScalarField>(alphaName_)
     && mesh.foundObject<volVectorField>(UName_)
    )
    {
        // Get alpha field
        const volScalarField& alpha1 =
            mesh.lookupObject<volScalarField>(alphaName_);

        // Mark all cells with alpha1 > alphaThreshold with monotonically
        // increasing counter
        labelList alphaCellIndexing(mesh.nCells(), -1);

        // Create a dynamic list containing alpha cells. Assume all cells will
        // be alpha cells to prevent any resizing
        dynamicLabelList alphaCells(mesh.nCells());

        const scalarField& alphaIn = alpha1.internalField();
        label nAlphaCells = 0;
        forAll (alphaIn, cellI)
        {
            if (alphaIn[cellI] > alphaThreshold_)
            {
                // This is alpha cell: mark it and append it to the list
                alphaCells.append(cellI);
                alphaCellIndexing[cellI] = nAlphaCells++; // Note postfix ++
            }
        }

        // Note on the algorithm: alphaCells are collected in increasing order
        // of cell indices, which means that we can loop through the cells in
        // the list, and for each cell, we can set face neighbours with higher
        // alpha cell index to the lower value of the owner. This means that all
        // alpha cells will be divided into uniquely labeled clusters after the
        // loop.

        // Get cell cells from the mesh
        const labelListList& meshCellCells = mesh.cellCells();

        // Loop control
        label nChangedAlphaIndices = 0;

        // Loop through all collected cells
        do
        {
            // Reset number of changed alpha indices and repeat until we have
            // identified all the clusters
            nChangedAlphaIndices = 0;

            forAll (alphaCells, i)
            {
                // Get cell index and alpha index
                const label& cellI = alphaCells[i];
                label& alphaCellI = alphaCellIndexing[cellI];

                // Loop through face neighbours
                const labelList& nbrCells = meshCellCells[cellI];
                forAll (nbrCells, nbrI)
                {
                    // Get cell index and alpha index
                    const label& nbrCellI = nbrCells[nbrI];
                    label& alphaCellNbr = alphaCellIndexing[nbrCellI];

                    if (alphaCellNbr != -1)
                    {
                        // This is neighbouring alpha cell, proceed

                        if (alphaCellNbr > alphaCellI)
                        {
                            // Alpha cell index of the neighbour is larger than
                            // alpha cell index of this cell, set the neighbour
                            // alpha index so that they are the same
                            alphaCellNbr = alphaCellI;

                            // Increment number of changed indices
                            ++nChangedAlphaIndices;
                        }
                        else if (alphaCellI > alphaCellNbr)
                        {
                            // Alpha cell index of the owner is larger than
                            // alpha cell index of this cell, set this cell's
                            // alpha index so that they are the same
                            alphaCellI = alphaCellNbr;

                            // Increment number of changed indices
                            ++nChangedAlphaIndices;
                        }
                        // else the alpha indices are already the same

                    } // End if the neighbouring cell is marked as alpha cell
                } // End loop for all neighbour cells
            } // End for all alpha cells

        } while (nChangedAlphaIndices != 0);

        // Loop through all alpha cells and collect unique indices
        label nUniqueClusters = 0;

        // Memory management
        {
            labelHashSet uniqueAlphaIDs(alphaCells.size());
            forAll (alphaCells, i)
            {
                // Get alpha index
                const label& alphaCellI = alphaCellIndexing[alphaCells[i]];

                // Insert into hash set
                uniqueAlphaIDs.insert(alphaCellI);
            }

            nUniqueClusters = uniqueAlphaIDs.size();
        }

        // We have ended up with nAlphaCells unique clusters of cells.
        // Create the list that will contain all the clusters
        PtrList<dynamicLabelList> alphaClusters(nUniqueClusters);

        // Store mapping from alphaCellIndexing to alphaCluster list:
        // key = alpha cell index
        // value = index into the PtrList
        Map<label> alphaCellIndexToClusterIndex(nUniqueClusters);

        // Set dummy value which indexes dummy alphaCellIndex value of -1 to
        // dummy cluster index of -1. Will be useful for generic treatment of
        // coupled boundaries below
        alphaCellIndexToClusterIndex.insert(-1, -1);

        // Helper variable for setting list members
        label listIndex = 0;

        // Loop through all alpha cells again
        forAll (alphaCells, i)
        {
            // Get cell index and alpha index
            const label& cellI = alphaCells[i];
            const label& alphaCellI = alphaCellIndexing[cellI];

            // Try inserting the new entry into the map
            if (alphaCellIndexToClusterIndex.insert(alphaCellI, listIndex))
            {
                // New entry inserted into the map

                // Set the cluster (dynamic list) at location listIndex
                // Note: A quite conservative (i.e. big) size estimate since we
                // can end up collecting a lot of cells in a cluster
                alphaClusters.set
                (
                    listIndex,
                    new dynamicLabelList
                    (
                        mesh.nCells()/16
                    )
                );

                // Append the cell to the cluster. Note: postfix increment to
                // increment the list index
                alphaClusters[listIndex++].append(cellI);
            }
            else
            {
                // The cluster with alphaCellI already exists, find the
                // corresponding list index

                // Get the mapping
                const label existingListIndex =
                    alphaCellIndexToClusterIndex[alphaCellI];

                // Append the cell to the existing cluster
                alphaClusters[existingListIndex].append(cellI);
            }
        }


        // At this point, we have all the clusters of connected alpha cells on
        // each processor separately. Let's calculate the volume (mass), centre
        // of volume (position), and velocity of each cluster
        List<clusterData> alphaClustersData(alphaClusters.size());

        // Get the velocity field and mesh data
        const vectorField& UIn =
            mesh.lookupObject<volVectorField>(UName_).internalField();

        const scalarField& VIn = mesh.V().field();
        const vectorField& CIn = mesh.C().internalField();

        forAll (alphaClusters, cI)
        {
            // Get the list of cells in this cluster
            const dynamicLabelList& clusterCells = alphaClusters[cI];

            // Loop through all the cells in this cluster and calculate the data
            scalar vc = 0;
            point  xc = point::zero;
            vector Uc = vector::zero;

            forAll (clusterCells, i)
            {
                // Get the cell and all the data
                const label& cellI = clusterCells[i];
                const scalar wetVolume = alphaIn[cellI]*VIn[cellI];
                const point& curC = CIn[cellI];
                const vector& curU = UIn[cellI];

                // Accumulate the data
                vc += wetVolume;
                xc += curC*wetVolume; // Note: do we need it more accurate?
                Uc += curU*wetVolume;
            }

            // Calculate the final position and velocity
            xc /= vc;
            Uc /= vc;

            // Set the data for this cluster
            alphaClustersData[cI].setData(vc, xc, Uc);
        }


        // Now the remaining part: we need to combine the data across coupled
        // boundaries such that the clusters touching each other on different
        // sides are merged together. All clusters near the boundary will be
        // shipped to master along with the connectivity and master will do all
        // the work. This part really needs to be "serialized" because certain
        // clusters can span across multiple processors and they should be
        // agglomerated one-by-one.
        // Algorithm:
        // 1. For all boundary faces, set local cluster indices from owner
        // 2. Swap the boundary face list to get the data from the other side
        // 3. Collect the following information:
        //    a) Create a global mapping that tells us which cluster on which
        //       processor speaks to which cluster on a different processor,
        //    b) Collect mapping from processor local cluster index to boundary
        //       cluster index,
        //    c) Collect cells that are on boundary into another list and
        //       "deactivate" the cluster in the original list to avoid
        //       duplication.
        // 4. Ship all the boundary droplets on all processors to master, along
        //    with the connectivity information.
        // 5. Master does data agglomeration.

        // Notes:

        // 1. If the parallel performance proves to be a bottleneck, this
        //    needs to be reorganized. Note that this won't be straightforward
        //    parallelisation using e.g. mapDistribute because a single cluster
        //    (i.e. droplet) can stretch across more than 2 processors and the
        //    agglomeration would need to happen by combining data from procI
        //    and procJ and then combining that to procK and so on.
        // 2. We still need to do this if we are running in serial because of
        //    cyclic patches (although it could be much simpler). I could do a
        //    slightly simplified procedure if (!Pstream::parRun()) but does
        //    anyone really care about minor loss of efficiency in serial?


        // PART 0: Data preparation

        // Create the cluster connectivity list. For my processor (first index),
        // we will collect a hash set containing unique pairs of cluster indices
        // across coupled patches, found on the other processor (second
        // index). The first entry in the pair will be a cluster index on my
        // (local) processor and the second entry will be a cluster index on
        // other (neighbouring) processor
        List<List<labelPairHashSet> > clusterConnectivity(Pstream::nProcs());
        forAll (clusterConnectivity, procI)
        {
            clusterConnectivity[procI].setSize(Pstream::nProcs());
        }
        List<labelPairHashSet>& myClusterConnectivity =
            clusterConnectivity[Pstream::myProcNo()]; // My proc part

        // Collect clusters that need parallel treatment (spanning across
        // multiple coupled boundaries)
        List<DynamicList<clusterData> > boundaryClusters(Pstream::nProcs());
        DynamicList<clusterData>& myBoundaryClusters =
            boundaryClusters[Pstream::myProcNo()]; // My proc part
        myBoundaryClusters.setCapacity(alphaClusters.size()/5); // Guessed size

        // Collect mapping from local cluster index to boundary clusters that
        // will be treated on master
        List<Map<label> > localToBoundary(Pstream::nProcs());
        Map<label>& myLocalToBoundary = localToBoundary[Pstream::myProcNo()];


        // PART 1: Set neighbouring cluster index on my own (owner) side

        // Communicate connectivity across coupled patches
        const label nFaces = mesh.nFaces();
        const label nInternalFaces = mesh.nInternalFaces();
        labelList neiClusterIndex(nFaces - nInternalFaces, -1);

        // Get face owners from polyMesh
        const labelList& meshFaceOwner = mesh.faceOwner();

        // Loop through boundary faces
        forAll (neiClusterIndex, i)
        {
            // Get face index and patch cell index
            const label faceI = i + nInternalFaces;
            const label& cellI = meshFaceOwner[faceI];

            // Get alpha cell index for this cell
            const label& acI = alphaCellIndexing[cellI];

            // Get corresponding cluster index
            const label& cI = alphaCellIndexToClusterIndex[acI];

            // Finally set the face value to the cluster index
            neiClusterIndex[i] = cI;
        }


        // PART 2: Swap the neighbouring cluster index across coupled patches

        // At this point, we have corresponding cluster indices for each face at
        // the boundary. Note that the cluster index is -1 where we don't
        // actually have clusters (droplets)

        // Swap the list
        syncTools::swapBoundaryFaceList(mesh, neiClusterIndex, false);


        // PART 3: Collect boundary clusters and necessary mappings

        // Get boundary mesh and alpha at the boundary
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
        const volScalarField::GeometricBoundaryField& alphab =
            alpha1.boundaryField();

        // Count number of global clusters on this processor
        label nBoundaryClusters = 0;

        // Now, we will compare the values from the two sides and fill in the
        // connectivity information and collect boundary clusters
        forAll (neiClusterIndex, i)
        {
            // Get face index and patch cell index
            const label faceI = i + nInternalFaces;
            const label& cellI = meshFaceOwner[faceI];

            // Get alpha cell index for this cell
            const label& acI = alphaCellIndexing[cellI];

            // Get corresponding cluster indices on my side and the other side
            const label& myCI = alphaCellIndexToClusterIndex[acI];
            const label& otherCI = neiClusterIndex[i];

            if ((myCI > -1) && (otherCI > -1))
            {
                // We have clusters on both sides of this face and these two
                // clusters have to be agglomerated

                // Get the patch index
                const label patchI = bMesh.whichPatch(faceI);

                // Get the fvPatchField of alpha in order to check whether the
                // fvPatchField is coupled and not the polyPatch
                const fvPatchScalarField& alphap = alphab[patchI];

                if (alphap.coupled())
                {
                    // This is a coupled patch that requires further treatment

                    if (isA<processorFvPatchScalarField>(alphap))
                    {
                        // Processor patch, cast
                        const processorFvPatchScalarField& procAlphap =
                            refCast<const processorFvPatchScalarField>(alphap);

                        // Get neighbouring processor index
                        const label nbrProcI = procAlphap.neighbProcNo();

                        // Part 3.a: Insert the connectivity my cluster has with
                        // the one on the other side
                        if
                        (
                            myClusterConnectivity[nbrProcI].insert
                            (
                                labelPair
                                (
                                    myCI,   // 1.: cluster ID on my side
                                    otherCI // 2.: cluster ID on other side
                                )
                            )
                        )
                        {
                            // PART 3.b: Unique pairs found and inserted. Set
                            // local to boundary mapping for my cluster index
                            if (myLocalToBoundary.insert(myCI, nBoundaryClusters))
                            {
                                // PART 3.c: New entry inserted into the map,
                                // append to the boundary list
                                myBoundaryClusters.append(alphaClustersData[myCI]);

                                // Deactivate the cluster in the local list
                                alphaClustersData[myCI].deactivate();

                                // Increment the counter
                                ++nBoundaryClusters;
                            }
                        }
                    }
                    else
                    {
                        // PART 3.a: This is an ordinary coupled patch, assume
                        // that both clusters are on the same processor and
                        // insert
                        if
                        (
                            myClusterConnectivity[Pstream::myProcNo()].insert
                            (
                                labelPair
                                (
                                    myCI,   // 1.: cluster ID (same proc)
                                    otherCI // 2.: other cluster ID (same proc)
                                )
                            )
                        )
                        {
                            // PART 3.b: Unique pairs found and inserted. Set
                            // local to boundary mapping for my cluster index
                            if (myLocalToBoundary.insert(myCI, nBoundaryClusters))
                            {
                                // PART 3.c: New entry inserted into the map,
                                // append to the boundary list
                                myBoundaryClusters.append(alphaClustersData[myCI]);

                                // Deactivate the cluster in the local list
                                alphaClustersData[myCI].deactivate();

                                // Increment the counter
                                ++nBoundaryClusters;
                            }
                        }

                    } // End if processor or ordinary coupled patch
                } // End if coupled patch
            } // End check whether there are clusters across faces
        } // End for all boundary faces


        // PART 4: We're all set with the connectivity information. Let's
        // communicate the connectivity, the clusters and the mapping to master
        Pstream::gatherList(clusterConnectivity);
        Pstream::gatherList(localToBoundary);
        Pstream::gatherList(boundaryClusters);


        // PART 5: Master does the agglomeration of clusters

        if (Pstream::master())
        {
            // Loop through all processors
            forAll (clusterConnectivity, procI)
            {
                // Get connectivity for this processor
                const List<labelPairHashSet>& procIConnectivity =
                    clusterConnectivity[procI];

                // Loop through remaining processors (and mine as well in case
                // we have cyclic patches)
                for (label procJ = procI; procJ < Pstream::nProcs(); ++procJ)
                {
                    // Get the connectivity from procI to procJ
                    const labelPairHashSet& procIprocJConnectivity =
                        procIConnectivity[procJ];

                    // Loop through the connectivity
                    forAllConstIter
                    (
                        labelPairHashSet,
                        procIprocJConnectivity,
                        iter
                    )
                    {
                        // Get cluster pairs
                        const labelPair& clusterIndices = iter.key();

                        // Get procI and procJ cluster indices
                        const label& cIProcI = clusterIndices.first();
                        const label& cIProcJ = clusterIndices.second();

                        // Get indices into boundary lists. Note: connectivity
                        // stored in procI only
                        const label& bcI = localToBoundary[procI][cIProcI];
                        const label& bcJ = localToBoundary[procJ][cIProcJ];

                        // Grab the clusters that need to be agglomerated
                        clusterData& boundaryClusterI =
                            boundaryClusters[procI][bcI];
                        clusterData& boundaryClusterJ =
                            boundaryClusters[procJ][bcJ];

                        // Calculate the agglomerated data
                        const scalar vc =
                            boundaryClusterI.volume()
                          + boundaryClusterJ.volume();
                        const point xc =
                        (
                            boundaryClusterI.centre()*boundaryClusterI.volume()
                          + boundaryClusterJ.centre()*boundaryClusterJ.volume()
                        )/vc;
                        const vector Uc =
                        (
                            boundaryClusterI.velocity()*
                            boundaryClusterI.volume()
                          + boundaryClusterJ.velocity()*
                            boundaryClusterJ.volume()
                        )/vc;

                        // Agglomerate clusterJ into clusterI
                        boundaryClusterI.setData(vc, xc, Uc);

                        // Deactivate clusterJ
                        boundaryClusterJ.deactivate();

                    } // End loop through all cross-processor connectivity
                } // End loop through remaining processors
            } // End loop through all processors
        } // End if master processor


        // Prepare for writing: count number of active clusters
        label nActiveClusters = 0;
        forAll (alphaClustersData, i)
        {
            if (alphaClustersData[i].active())
            {
                ++nActiveClusters;
            }
        }

        // Count boundary clusters only for master since master handles them
        if (Pstream::master())
        {
            forAll (boundaryClusters, procI)
            {
                const DynamicList<clusterData>& bClusters =
                    boundaryClusters[procI];

                forAll (bClusters, i)
                {
                    if (bClusters[i].active())
                    {
                        ++nActiveClusters;
                    }
                }
            }
        }

        // Now we have all the data we need. Let's dump it into files. Note:
        // each processor writes separate file with the following format:
        // caseDir/dropletStatistics/t0.01-proc0.dat
        if (nActiveClusters > 0)
        {
            // Get file name
            fileName ofName;

            if (Pstream::parRun())
            {
                ofName =
                    mesh.time().rootPath()/
                    mesh.time().caseName()/
                    ".."/
                    dirName_/
                    (
                        "t"
                      + mesh.time().timeName()
                      + "-proc"
                      + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    );
            }
            else
            {
                ofName =
                    mesh.time().rootPath()/
                    mesh.time().caseName()/
                    dirName_/
                    (
                        "t"
                      + mesh.time().timeName()
                      + "-proc"
                      + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    );
            }

            // Create stream for writing
            OFstream ofs(ofName);

            // Write header
            ofs << "# Position, velocity, diameter" << endl;

            // Loop through all local data and write to file
            forAll (alphaClustersData, i)
            {
                // Get current cluster
                const clusterData& cd = alphaClustersData[i];

                if (cd.active())
                {
                    // Check whether the diameter is within the bounds
                    if
                    (
                        cd.diameterWithinBounds(minDiameter_, maxDiameter_)
                    )
                    {
                        ofs << cd.centre() << tab
                            << cd.velocity() << tab
                            << cd.effectiveDiameter() << endl;
                    }
                }
            }

            // Master has additional job to do: write all boundary
            // droplets. Note that these boundary droplets are all "identified"
            // by master processor. This doesn't mean that these droplets are
            // actually on the master processor.
            if (Pstream::master())
            {
                // Loop through all processors
                forAll (boundaryClusters, procI)
                {
                    // Get list of boundary clusters owned by this processor
                    const DynamicList<clusterData>& bClusters =
                        boundaryClusters[procI];

                    // Loop through all clusters
                    forAll (bClusters, i)
                    {
                        // Get this boundary cluster
                        const clusterData& bcd = bClusters[i];

                        if (bcd.active())
                        {
                            // Check whether the diameter is within the bounds
                            if
                            (
                                bcd.diameterWithinBounds
                                (
                                    minDiameter_,
                                    maxDiameter_
                                )
                            )
                            {
                                ofs << bcd.centre() << tab
                                    << bcd.velocity() << tab
                                    << bcd.effectiveDiameter() << nl;
                            }
                        }
                    } // End loop over all boundary clusters
                } // End loop over all processors
            } // End if master processor

            // Flush the stream
            ofs.flush();

        } // End if there are active clusters
    }
    else
    {
        InfoIn("bool dropletStatistics::execute()")
            << "Did not find volume fraction field: " << alphaName_
            << " and/or velocity field: " << UName_
            << nl << "Doing nothing..."
            << endl;
    }
}


// ************************************************************************* //
