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
    Foam::deflatedPCG

Description
    Preconditioned conjugate gradient solver for symmetric lduMatrices
    using a run-time selectable preconditioner.

    Matrix deflation wass added. This solver only works for singular matrices
    like derived from Laplacian with only Neumann-BC and WITHOUT reference set.

SourceFiles
    deflatedPCG.C

\*---------------------------------------------------------------------------*/

#ifndef deflatedPCG_H
#define deflatedPCG_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deflatedPCG Declaration
\*---------------------------------------------------------------------------*/

class deflatedPCG
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        deflatedPCG(const deflatedPCG&);

        //- Disallow default bitwise assignment
        void operator=(const deflatedPCG&);


public:

    //- Runtime type information
    TypeName("deflatedPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        deflatedPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );


    // Destructor

        virtual ~deflatedPCG()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual lduSolverPerformance solve
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
