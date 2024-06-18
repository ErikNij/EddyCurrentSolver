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

Class
    Foam::vectorMetisDecomp

Description
    Metis domain decomposition with vector weights

SourceFiles
    vectorMetisDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef vectorMetisDecomp_H
#define vectorMetisDecomp_H

#include "decompositionMethod.H"

extern "C"
{
#   include "metis.h"
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class vectorMetisDecomp Declaration
\*---------------------------------------------------------------------------*/

class vectorMetisDecomp
:
    public decompositionMethod
{
    // Private structs

        struct PartGraphData
        {
            const int nvtxs;
            const int nadjs;
            const List<int>& adjncy;
            const List<int>& xadj;
            const int nparts;

            mutable int ncon;
            mutable List<int> vwgt;
            mutable List<int> adjwgt;
            mutable List<real_t> tpwgts;
            mutable List<real_t> ubvec;

            PartGraphData
            (
                const List<int>& adjncy_,
                const List<int>& xadj_,
                int nparts_
            )
            :
                nvtxs(xadj_.size() - 1),
                nadjs(adjncy_.size()/2),
                adjncy(adjncy_),
                xadj(xadj_),
                nparts(nparts_),
                ncon(1)
            {}
        };


    // Private data

        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        vectorMetisDecomp(const vectorMetisDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const vectorMetisDecomp&);

        //- Calculate Metis graph data for partitioning.
        //  Various settings are processed here
        PartGraphData calcPartGraphData
        (
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& weights,
            const labelList& agglom = labelList()
        );

        //- Basic checks for graph data
        bool checkPartGraphData
        (
            const word& invokedFrom,
            const PartGraphData& pgd
        );

        //- Partiton graph
        int decompose
        (
            const PartGraphData& pgd,
            List<int>& part
        );


public:

    //- Runtime type information
    TypeName("vectorMetis");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        vectorMetisDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~vectorMetisDecomp()
        {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // Metis does not know about proc boundaries
            return false;
        }

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField& regionPoints,
            const scalarField& regionWeights
        );

        //- Same but with uniform weights
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField& regionPoints
        )
        {
            return decompose
            (
                agglom,
                regionPoints,
                scalarField(regionPoints.size(), 1.0)
            );
        }

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
