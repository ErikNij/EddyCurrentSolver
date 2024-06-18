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
    Foam::newLiquidFilmFvMotionSolver

Description
    Mesh motion solver for an fvMesh.

SourceFiles
    newLiquidFilmFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef newLiquidFilmFvMotionSolver_H
#define newLiquidFilmFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class newLiquidFilmFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class newLiquidFilmFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- Point motion field
        mutable pointVectorField pointMotionU_;

        //- Point to point addressing
        labelList pointPoint_;

        //- Top plane y-coordinate
        scalar max_;

        //- Bottom plane y-coordinate
        scalar min_;

        //- Motion direction
        vector motionDirection_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        newLiquidFilmFvMotionSolver
        (
            const newLiquidFilmFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const newLiquidFilmFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("newLiquidFilm");


    // Constructors

        //- Construct from polyMesh and data stream
        newLiquidFilmFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        virtual ~newLiquidFilmFvMotionSolver();


    // Member Functions

        //- Return reference to the point motion velocity field
        pointVectorField& pointMotionU()
        {
            return pointMotionU_;
        }

        //- Return const reference to the point motion velocity field
        const pointVectorField& pointMotionU() const
        {
            return pointMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
