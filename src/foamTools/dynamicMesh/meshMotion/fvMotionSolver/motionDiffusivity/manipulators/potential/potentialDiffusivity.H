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
    Foam::potentialDiffusivity

Description
    Mesh motion diffusivity manipulator which returns the TODO
    of the given diffusivity D.

SourceFiles
    potentialDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef potentialDiffusivity_H
#define potentialDiffusivity_H

#include "motionDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class potentialDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class potentialDiffusivity
:
    public motionDiffusivity
{

    // Private data

        //- Exponent
        scalar alpha_;

        //- Basic diffusivity
        autoPtr<motionDiffusivity> basicDiffusivityPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        potentialDiffusivity(const potentialDiffusivity&);

        //- Disallow default bitwise assignment
        void operator=(const potentialDiffusivity&);


public:

    //- Runtime type information
    TypeName("potential");


    // Constructors

        //- Construct for the given fvMotionSolver and data Istream
        potentialDiffusivity
        (
            const fvMotionSolver& mSolver,
            Istream& mdData
        );


    // Destructor

        virtual ~potentialDiffusivity();


    // Member Functions

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const;

        //- Correct the motion diffusivity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
