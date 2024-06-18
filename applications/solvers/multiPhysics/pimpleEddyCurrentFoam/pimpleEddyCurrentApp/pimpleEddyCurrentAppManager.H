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
    Foam::pimpleEddyCurrentApp::Manager

Description
    ...

SourceFiles
    pimpleEddyCurrentAppManager.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef pimpleEddyCurrentAppManager_H
#define pimpleEddyCurrentAppManager_H

#include "pimpleEddyCurrentAppInclude.H"
#include "pimpleEddyCurrentAppControl.H"

#include "pimpleAppManager.H"
#include "eddyCurrentAppManager.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace pimpleEddyCurrentApp
{

/*---------------------------------------------------------------------------*\
               Class pimpleEddyCurrentApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public regionFvSolverManager
{
public:

    // Public typedefs and forward declarations

        typedef pimpleEddyCurrentApp::Control Control;
        class Settings;
        class Storage;

    // Solver Manager Construction

        SM_GLOBALCONTROL_ADD_P1
        (
            Control,
            mesh()
        );

        SM_GLOBALSETTINGS_NEW();
        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();
        SM_GLOBALSTORAGE_ADD();


        SM_REGIONS_NEW();
        SM_REGIONS_ADD();


private:

    // Private data

        //- Pointer for manager of pimpleApp
        mutable autoPtr<pimpleApp::Manager> pimpleAppManager_;

        //- Pointer for manager of eddyCurrentApp
        mutable autoPtr<eddyCurrentApp::Manager> eddyCurrentAppManager_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);


public:

    // Runtime type information
    TypeName("pimpleEddyCurrentApp::Manager");

    // Constructor

        Manager
        (
            const argList& args,
            Time& time,
            regionFvMesh& mesh,
            bool master = true
        );

    // Destructor

        virtual ~Manager();

    // Public Member Functions

        // Managers

            //- pimpleApp
            pimpleApp::Manager& pimpleAppManager() const;

            //- eddyCurrentApp
            eddyCurrentApp::Manager& eddyCurrentAppManager() const;

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

            //- Read sequence
            virtual void read() const;

            //- Init sequence
            virtual void init() const;

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pimpleEddyCurrentApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pimpleEddyCurrentAppManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
