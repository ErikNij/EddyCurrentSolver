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

\*---------------------------------------------------------------------------*/

#include "vectorMetisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vectorMetisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        vectorMetisDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vectorMetisDecomp::PartGraphData Foam::vectorMetisDecomp::calcPartGraphData
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& weights,
    const labelList& agglom
)
{
    PartGraphData pgd(adjncy, xadj, nProcessors_);

    // Check for externally provided weights and if so initialise vwgt
    scalar minWeights = gMin(weights);
    if (weights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "vectorMetisDecomp::calcPartGraphData"
                "(const List<int>&, const List<int>&,"
                " const scalarField&, const labelList&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        // Convert to integers.
        pgd.vwgt.setSize(weights.size());
        forAll(pgd.vwgt, i)
        {
            pgd.vwgt[i] = int(weights[i]/minWeights);
        }
    }

    // Check for user supplied weights
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");

        // Number of cells
        int nCells = (agglom.size() > 0) ? agglom.size() : pgd.nvtxs;

        word weightsFile;
        if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "vectorMetisDecomp : Using cell-based weights" << endl;

            IOList<int> cellIOWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );

            if (cellIOWeights.size() != nCells)
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::calcPartGraphData"
                    "(const List<int>&, const List<int>&,"
                    " const scalarField&, const labelList&)"
                )   << "Number of cell weights " << cellIOWeights.size()
                    << " does not equal number of cells " << nCells
                    << exit(FatalError);
            }

            if (agglom.size() > 0)
            {
                // Initialise vertice weights if necessary
                if (pgd.vwgt.size() == 0)
                {
                    pgd.vwgt.setSize(pgd.nvtxs);
                }

                forAll(cellIOWeights, celli)
                {
                    pgd.vwgt[agglom[celli]] = cellIOWeights[celli];
                }
            }
            else
            {
                pgd.vwgt.transfer(cellIOWeights);
            }
        }

        if (metisCoeffs.lookupOrDefault("regionBalancing", false))
        {
            Info<< "vectorMetisDecomp : Using region-balancing"
                << " multi-constraint feature" << endl;

            wordList regionCellSetNames =
                metisCoeffs.lookupOrDefault("regionCellSets", wordList());

            if (regionCellSetNames.size() == 0)
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::calcPartGraphData"
                    "(const List<int>&, const List<int>&,"
                    " const scalarField&, const labelList&)"
                )   << "Feature region-balancing is used but list of region"
                    << " cell set names 'regionCellSets' is empty."
                    << exit(FatalError);
            }

            Info<< "vectorMetisDecomp : Balancing regions from cell sets"
                << regionCellSetNames << endl;

            Field<real_t> regionWeightFactors(regionCellSetNames.size(), 1.0);
            if (metisCoeffs.readIfPresent("regionWeightFactors", regionWeightFactors))
            {
                if (regionWeightFactors.size() != regionCellSetNames.size())
                {
                    // Use one value for all scale factors
                    if (regionWeightFactors.size() == 1)
                    {
                        real_t singleScale = regionWeightFactors[0];

                        regionWeightFactors.setSize
                        (
                            regionCellSetNames.size(),
                            singleScale
                        );
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "vectorMetisDecomp::calcPartGraphData"
                            "(const List<int>&, const List<int>&,"
                            " const scalarField&, const labelList&)"
                        )   << "Number of region weight scale factors "
                            << regionWeightFactors.size() << " is not 1 or"
                            << " does not equal number of given region cell sets "
                            << regionCellSetNames.size()
                            << exit(FatalError);
                    }
                }

                Info<< "vectorMetisDecomp : Scaling region weights according to "
                    << "factors " << regionWeightFactors << endl;
            }

            // Each region gets its own multi-constraint vector index
            pgd.ncon += regionCellSetNames.size();

            // Cache current weights as global or create default ones
            List<int> globalWeights(pgd.vwgt);
            if (globalWeights.size() == 0)
            {
                globalWeights.setSize(pgd.nvtxs, 1);
            }

            // Resize weights to ncon*nvtxs and initialise with zero
            pgd.vwgt.setSize(pgd.ncon*pgd.nvtxs, 0);

            // Insert global weights as first multi-constraint vector index
            forAll (globalWeights, i)
            {
                int ci = pgd.ncon*i;

                pgd.vwgt[ci] = globalWeights[i];
            }

            // Insert weights on additional multi-constraint vector indices
            // for each named region
            forAll (regionCellSetNames, ni)
            {
                cellSet curCellSet
                (
                    mesh_,
                    regionCellSetNames[ni],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                List<int> cellLabels = curCellSet.toc();

                int scale(regionWeightFactors[ni]*pow(1.0*nCells/cellLabels.size(),2));

                if (agglom.size() > 0)
                {
                    forAll(cellLabels, celli)
                    {
                        int i = agglom[cellLabels[celli]];

                        int cI = pgd.ncon*i + ni + 1;

                        pgd.vwgt[cI] = scale*globalWeights[i];
                    }
                }
                else
                {
                    forAll(cellLabels, celli)
                    {
                        int i = cellLabels[celli];

                        int cI = pgd.ncon*i + ni + 1;

                        pgd.vwgt[cI] = scale*globalWeights[i];
                    }
                }
            }
        }

        Field<real_t> processorWeights;
        if (metisCoeffs.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            if (processorWeights.size() != pgd.nparts)
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::calcPartGraphData"
                    "(const List<int>&, const List<int>&,"
                    " const scalarField&, const labelList&)"
                )   << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number graph partitions" << pgd.nparts
                    << exit(FatalError);
            }

            Info<< "vectorMetisDecomp : Using processor weights "
                << processorWeights << endl;

            // Equal partiton weights for each constraint
            pgd.tpwgts.setSize(pgd.nparts*pgd.ncon);
            for (int i = 0; i < pgd.nparts; i++)
            {
                for (int j = 0; j < pgd.ncon; j++)
                {
                    pgd.tpwgts[pgd.ncon*i + j] = processorWeights[i];
                }
            }
        }

        Field<real_t> imbTolerances;
        if (metisCoeffs.readIfPresent("imbalanceTolerances", imbTolerances))
        {
            if (imbTolerances.size() != pgd.ncon)
            {
                // Use one value for all tolerances
                if (imbTolerances.size() == 1)
                {
                    real_t singleTol = imbTolerances[0];
                    imbTolerances.setSize(pgd.ncon, singleTol);
                }
                else
                {
                    FatalErrorIn
                    (
                        "vectorMetisDecomp::calcPartGraphData"
                        "(const List<int>&, const List<int>&,"
                        " const scalarField&, const labelList&)"
                    )   << "Number of graph load imbalance tolerances "
                        << imbTolerances.size() << " is not 1 or"
                        << " does not equal number of graph constraints "
                        << pgd.ncon
                        << exit(FatalError);
                }
            }

            if (min(imbTolerances) < 1.0)
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::calcPartGraphData"
                    "(const List<int>&, const List<int>&,"
                    " const scalarField&, const labelList&)"
                )   << "Load imbalance tolerances must be at least 1.0"
                    << exit(FatalError);
            }

            Info<< "vectorMetisDecomp : Using load imbalance tolerances "
                << imbTolerances << endl;

            pgd.ubvec.transfer(imbTolerances);
        }
    }

    return pgd;
}


bool Foam::vectorMetisDecomp::checkPartGraphData
(
    const word& name,
    const PartGraphData& pgd
)
{
    // Sanity check for vertice weights
    if ((pgd.vwgt.size() > 0) && (pgd.vwgt.size() != pgd.ncon*pgd.nvtxs))
    {
        FatalErrorIn(name)
            << "Number of graph vertice weights " << pgd.vwgt.size()
            << " does not equal number of graph constraints " << pgd.ncon
            << " times number of graph vertices " << pgd.nvtxs
            << exit(FatalError);
    }

    // Sanity check for edge weights
    if ((pgd.adjwgt.size() > 0) && (pgd.adjwgt.size() != pgd.nadjs))
    {
        FatalErrorIn(name)
            << "Number of graph edge weights " << pgd.adjwgt.size()
            << " does not equal number of graph edges " << pgd.nadjs
            << exit(FatalError);
    }

    // Sanity check for partition weights
    if ((pgd.tpwgts.size() > 0) && (pgd.tpwgts.size() != pgd.ncon*pgd.nparts))
    {
        FatalErrorIn(name)
            << "Number of graph partition weights " << pgd.tpwgts.size()
            << " does not equal number of constraints " << pgd.ncon
            << " times number of graph partitions " << pgd.nparts
            << exit(FatalError);
    }

    // Sanity check for load imblance tolerances
    if ((pgd.ubvec.size() > 0) && (pgd.ubvec.size() != pgd.ncon))
    {
        FatalErrorIn(name)
            << "Number of graph load imbalance tolerances" << pgd.ubvec.size()
            << " does not equal number of graph constraints " << pgd.ncon
            << exit(FatalError);
    }

    return true;
}


int Foam::vectorMetisDecomp::decompose
(
    const PartGraphData& pgd,
    List<int>& part
)
{
    // Sanity checks for graph data
    checkPartGraphData
    (
        "vectorMetisDecomp::decompose(const PartGraphData&, List<int>&)",
        pgd
    );

    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("k-way");

    // decomposition options. 0 = use defaults
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");

        if (metisCoeffs.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::decompose"
                    "(const PartGraphData&, List<int>&)"
                )   << "Method " << method
                    << " in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "vectorMetisDecomp : Using Metis method " << method
                << nl << endl;
        }

        List<int> mOptions;
        if (metisCoeffs.readIfPresent("options", mOptions))
        {
            if (mOptions.size() != METIS_NOPTIONS)
            {
                FatalErrorIn
                (
                    "vectorMetisDecomp::decompose"
                    "(const PartGraphData&, List<int>&)"
                )   << "Number of options in metisCoeffs in dictionary "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            forAll(mOptions, i)
            {
                options[i] = mOptions[i];
            }

            Info<< "vectorMetisDecomp : Using Metis options " << mOptions
                << nl << endl;
        }
    }

    // nvtxs: num vertices in graph
    int nvtxs = pgd.nvtxs;

    // npart: num partitions
    int npart = pgd.nparts;

    // objval: number of cut edges
    int edgeCut = 0;

    // part: partition vector
    part.setSize(pgd.nvtxs);

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &nvtxs,             // nvtxs: num vertices in graph
            &pgd.ncon,          // ncon: num balancing constraints
            const_cast<List<int>&>(pgd.xadj).begin(), // xadj: indexing into adjncy
            const_cast<List<int>&>(pgd.adjncy).begin(), // adjncy: neighbour info
            pgd.vwgt.begin(),   // vwgt: vertex weights
            NULL,               // vsize: total communication vol
            pgd.adjwgt.begin(), // adjwgt: edge weights
            &npart,             // nparts: num partitions
            pgd.tpwgts.begin(), // tpwgts: partition weights
            pgd.ubvec.begin(),  // ubvec: partition imbalance
            options,
            &edgeCut,           // objval: number of cut edges
            part.begin()        // part: partition vector
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &nvtxs,             // nvtxs: num vertices in graph
            &pgd.ncon,          // ncon: num balancing constraints
            const_cast<List<int>&>(pgd.xadj).begin(), // xadj: indexing into adjncy
            const_cast<List<int>&>(pgd.adjncy).begin(), // adjncy: neighbour info
            pgd.vwgt.begin(),   // vwgt: vertex weights
            NULL,               // vsize: total communication vol
            pgd.adjwgt.begin(), // adjwgt: edge weights
            &npart,             // nparts: num partitions
            pgd.tpwgts.begin(), // tpwgts: partition weights
            pgd.ubvec.begin(),  // ubvec: partition imbalance
            options,
            &edgeCut,
            part.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vectorMetisDecomp::vectorMetisDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::vectorMetisDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "vectorMetisDecomp::decompose(const pointField&,const scalarField&)"
        )   << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    List<int> adjncy;
    List<int> xadj;
    calcCSR
    (
        mesh_,
        adjncy,
        xadj
    );

    // Decompose using weights
    List<int> finalDecomp;
    decompose
    (
        calcPartGraphData(adjncy, xadj, pointWeights),
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


Foam::labelList Foam::vectorMetisDecomp::decompose
(
    const labelList& agglom,
    const pointField& regionPoints,
    const scalarField& regionWeights
)
{
    if (agglom.size() != mesh_.nCells())
    {
        FatalErrorIn
        (
            "vectorMetisDecomp::decompose"
            "(const labelList&, const pointField&, const scalarField&)"
        )   << "Size of cell-to-coarse map " << agglom.size()
            << " differs from number of cells in mesh " << mesh_.nCells()
            << exit(FatalError);
    }

    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    List<int> adjncy;
    List<int> xadj;
    {
        // Get cellCells on coarse mesh.
        labelListList cellCells;
        calcCellCells
        (
            mesh_,
            agglom,
            regionPoints.size(),
            cellCells
        );

        calcCSR(cellCells, adjncy, xadj);
    }

    // Decompose using weights
    List<int> finalDecomp;
    decompose
    (
        calcPartGraphData(adjncy, xadj, regionWeights, agglom),
        finalDecomp
    );


    // Rework back into decomposition for original mesh_
    labelList decomp(agglom.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[agglom[i]];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


Foam::labelList Foam::vectorMetisDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorIn
        (
            "vectorMetisDecomp::decompose"
            "(const pointField&, const labelListList&, const scalarField&)"
        )   << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }


    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli

    List<int> adjncy;
    List<int> xadj;
    calcCSR(globalCellCells, adjncy, xadj);

    // Decompose using weights
    List<int> finalDecomp;
    decompose
    (
        calcPartGraphData(adjncy, xadj, cellWeights),
        finalDecomp
    );

    // Copy back to labelList
    labelList decomp(finalDecomp.size());

    forAll(decomp, i)
    {
        decomp[i] = finalDecomp[i];
    }

    fixCyclics(mesh_, decomp);

    return decomp;
}


// ************************************************************************* //
