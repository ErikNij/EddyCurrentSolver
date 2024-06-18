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
    Foam::refVelocityPseudoSolidFvMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    pseudo-solid motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable cellMotionU.
    Pseudo-solid motion equation is solved on inital mesh.

    Note: Diffusion field is used as Lame's constant \mu

SourceFiles
    refVelocityPseudoSolidFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef refVelocityPseudoSolidFvMotionSolver_H
#define refVelocityPseudoSolidFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
              Class refVelocityPseudoSolidFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class refVelocityPseudoSolidFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- Point motion field
        mutable pointVectorField pointMotionU_;

        //- Cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Use leastSquares volToPoint interpolaltion
        Switch leastSquaresVolPoint_;

        //- Poisson's ration
        scalar nu_;

        //- Number of correctors in solution procedure of pseudo-solid eqn
        label nCorrectors_;

        //- Convergence tolerance for solution of pseudo-solid eqn
        scalar convergenceTolerance_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        refVelocityPseudoSolidFvMotionSolver
        (
            const refVelocityPseudoSolidFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const refVelocityPseudoSolidFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("refVelocityPseudoSolid");


    // Constructors

        //- Construct from polyMesh and data stream
        refVelocityPseudoSolidFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        virtual ~refVelocityPseudoSolidFvMotionSolver();


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

        //- Return reference to the cell motion velocity field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell motion velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
