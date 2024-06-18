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
    Foam::solverPerformanceData

Description
    Class to hold SolverPerformance data for mesh database

SourceFiles
    solverPerformanceData.C

\*---------------------------------------------------------------------------*/

#ifndef solverPerformanceData_H
#define solverPerformanceData_H

#include "pTraits.H"
#include "IOstreams.H"
#include "scalarField.H"

#include "lduMatrix.H"
#include "BlockLduMatrix.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class solverPerformanceData;

Istream& operator>>
(
    Istream&,
    solverPerformanceData&
);

Ostream& operator<<
(
    Ostream&,
    const solverPerformanceData&
);

/*---------------------------------------------------------------------------*\
                    Class solverPerformanceData Declaration
\*---------------------------------------------------------------------------*/

class solverPerformanceData
{
protected:

    // Protected data

        //- Solver name
        word solverName_;

        //- Field name
        word fieldName_;

        //- Initial residual as field
        scalarField initialResidual_;

        //- Final residual as field
        scalarField finalResidual_;

        //- Number of iterations
        label nIterations_;

        //- Has the solver converged?
        bool converged_;

        //- Is the matrix singular?
        bool singular_;


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("solverPerformanceData");

    // Constructors

        //- Null constructor
        solverPerformanceData()
        :
            solverName_(),
            fieldName_(),
            initialResidual_(),
            finalResidual_(),
            nIterations_(-1),
            converged_(false),
            singular_(false)
        {}

        //- Construct with solver and field name
        solverPerformanceData
        (
            const word& solverName,
            const word& fieldName,
            const scalarField& iRes,
            const scalarField& fRes,
            const label nIter,
            const bool converged,
            const bool singular
        )
        :
            solverName_(solverName),
            fieldName_(fieldName),
            initialResidual_(iRes),
            finalResidual_(fRes),
            nIterations_(nIter),
            converged_(converged),
            singular_(singular)
        {}

        // Construct from BlockSolverPerformance<Type>
        template<class Type>
        explicit solverPerformanceData
        (
            const BlockSolverPerformance<Type>& sp
        )
        :
            solverName_(sp.solverName()),
            fieldName_(sp.fieldName()),
            initialResidual_(pTraits<Type>::nComponents, 0),
            finalResidual_(pTraits<Type>::nComponents, 0),
            nIterations_(sp.nIterations()),
            converged_(sp.converged()),
            singular_(sp.singular())
        {
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                initialResidual_[cmpt] = sp.initialResidual().component(cmpt);
                finalResidual_[cmpt] = sp.finalResidual().component(cmpt);
            }
        }

    // Member functions

        //- Return solver name
        inline const word& solverName() const
        {
            return solverName_;
        }

        //- Return solver name
        inline word& solverName()
        {
            return solverName_;
        }

        //- Return field name
        inline const word& fieldName() const
        {
            return fieldName_;
        }

        //- Return initial residual
        inline const scalarField& initialResidual() const
        {
            return initialResidual_;
        }

        //- Return initial residual
        inline scalarField& initialResidual()
        {
            return initialResidual_;
        }

        //- Return min initial residual
        inline scalar minInitialResidual() const
        {
            return min(initialResidual_);
        }

        //- Return max initial residual
        inline scalar maxInitialResidual() const
        {
            return max(initialResidual_);
        }

        //- Return final residual
        inline const scalarField& finalResidual() const
        {
            return finalResidual_;
        }

        //- Return final residual
        inline scalarField& finalResidual()
        {
            return finalResidual_;
        }

        //- Return min final residual
        inline scalar minFinalResidual() const
        {
            return min(finalResidual_);
        }

        //- Return max final residual
        inline scalar maxFinalResidual() const
        {
            return max(finalResidual_);
        }

        //- Return number of iterations
        inline label nIterations() const
        {
            return nIterations_;
        }

        //- Return number of iterations
        inline label& nIterations()
        {
            return nIterations_;
        }

        //- Has the solver converged?
        inline bool converged() const
        {
            return converged_;
        }

        //- Has the solver converged?
        inline bool& converged()
        {
            return converged_;
        }

        //- Is the matrix singular?
        inline bool singular() const
        {
            return singular_;
        }

        //- Is the matrix singular?
        inline bool& singular()
        {
            return singular_;
        }

    // Member Operators

        bool operator!=(const solverPerformanceData&) const;

    // IOstream Operators

        friend Istream& operator>>
        (
            Istream&,
            solverPerformanceData&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const solverPerformanceData&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
