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
    Foam::inverseDeltaDiffusivity

Description
    Inverse cell-volume motion diffusivity.

SourceFiles
    inverseDeltaDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef inverseDeltaDiffusivity_H
#define inverseDeltaDiffusivity_H

#include "uniformDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inverseDeltaDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class inverseDeltaDiffusivity
:
    public uniformDiffusivity
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        inverseDeltaDiffusivity(const inverseDeltaDiffusivity&);

        //- Disallow default bitwise assignment
        void operator=(const inverseDeltaDiffusivity&);


public:

    //- Runtime type information
    TypeName("inverseDelta");


    // Constructors

        //- Construct for the given fvMotionSolver and data Istream
        inverseDeltaDiffusivity
        (
            const fvMotionSolver& mSolver,
            Istream& mdData
        );


    // Destructor

        virtual ~inverseDeltaDiffusivity();


    // Member Functions

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
