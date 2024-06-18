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

Description
    ...

SourceFiles
    interTrackEddyCurrentAppControlI.H
    interTrackEddyCurrentAppControl.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef interTrackEddyCurrentAppControl_H
#define interTrackEddyCurrentAppControl_H

#include "interTrackEddyCurrentAppInclude.H"

#include "solutionControl.H"
#include "solverPerformanceData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interTrackEddyCurrentApp
{

/*---------------------------------------------------------------------------*\
                   Class interTrackEddyCurrentApp::Control Declaration
\*---------------------------------------------------------------------------*/

class Control
:
    public solutionControl
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        Control(const Control&);

        //- Disallow default bitwise assignment
        void operator=(const Control&);


protected:

    // Protected data

        //- Mesh
        const regionFvMesh& mesh_;

    // Protected Member Functions

        //- Return true if an electromagnetic update is not necessary
        virtual bool criteriaSatisfied();


public:

    //- Run-time type information
    TypeName("interTrackEddyCurrentApp::Control");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        Control
        (
            regionFvMesh& mesh,
            const word& dictName="MULTIPHYSICS"
        );


    //- Destructor
    virtual ~Control();


    // Member Functions

        // Evolution

            //- Main control loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interTrackEddyCurrentApp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interTrackEddyCurrentAppControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
