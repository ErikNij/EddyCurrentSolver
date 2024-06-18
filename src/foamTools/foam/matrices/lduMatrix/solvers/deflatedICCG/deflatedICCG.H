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
    Foam::deflatedICCG

Description
    Incomplete Cholesky preconditioned CG solver derived from the general
    preconditioned CG solver PCG but with the choice of preconditioner
    pre-selected.

    Matrix deflation wass added. This solver only works for singular matrices
    like derived from Laplacian with only Neumann-BC and WITHOUT reference set.

Deprecated
    This solver is present for backward-compatibility and the PCG solver
    should be used for preference.

SourceFiles
    deflatedICCG.C

\*---------------------------------------------------------------------------*/

#ifndef deflatedICCG_H
#define deflatedICCG_H

#include "deflatedPCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class deflatedICCG Declaration
\*---------------------------------------------------------------------------*/

class deflatedICCG
:
    public deflatedPCG
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        deflatedICCG(const deflatedICCG&);

        //- Disallow default bitwise assignment
        void operator=(const deflatedICCG&);

public:

        //- Return the dictionary constructed from the components.
        //  Needed for backward compatibility
        static dictionary solverDict
        (
            const scalar tol,
            const scalar relTol
        );

        //- Return the dictionary constructed from the old-style data-stream.
        //  Needed for backward compatibility
        static dictionary solverDict(Istream&);


    //- Runtime type information
    TypeName("deflatedICCG");


    // Constructors

        //- Construct from matrix components and solver data stream
        deflatedICCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );

        //- Construct from matrix components and tolerances
        deflatedICCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const scalar tolerance,
            const scalar relTol = 0.0
        );

    // Destructor

        virtual ~deflatedICCG()
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
