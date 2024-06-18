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
    Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager

Description
    ...

SourceFiles
    buoyantBoussinesqPimpleEddyCurrentAppManager.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef buoyantBoussinesqPimpleEddyCurrentAppManager_H
#define buoyantBoussinesqPimpleEddyCurrentAppManager_H

#include "buoyantBoussinesqPimpleEddyCurrentAppInclude.H"
#include "buoyantBoussinesqPimpleEddyCurrentAppControl.H"

#include "buoyantBoussinesqPimpleAppManager.H"
#include "eddyCurrentAppManager.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace buoyantBoussinesqPimpleEddyCurrentApp
{

/*---------------------------------------------------------------------------*\
        Class buoyantBoussinesqPimpleEddyCurrentApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public regionFvSolverManager
{
public:

    // Public typedefs and forward declarations

        typedef buoyantBoussinesqPimpleEddyCurrentApp::Control Control;
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

            //- Temperature
            SM_GLOBALITEM_ADD(regionVolScalarField, T);

            //- Thermal conductivity
            SM_GLOBALITEM_ADD(regionVolScalarField, lambda);

            //- Density times specific heat capacity
            SM_GLOBALITEM_ADD(regionVolScalarField, rhoCp);

// TODO: Substitute U with phi
            //- phi
//             SM_GLOBALITEM_ADD(regionSurfaceScalarField, phi);
            SM_GLOBALITEM_ADD(regionVolVectorField, U);

        SM_GLOBALSTORAGE_ADD();


        SM_REGIONS_NEW();
        SM_REGIONS_ADD();


private:

    // Private data

        //- Pointer for manager of buoyantBoussinesqPimpleApp
        mutable autoPtr<buoyantBoussinesqPimpleApp::Manager> buoyantBoussinesqPimpleAppManager_;

        //- Pointer for manager of eddyCurrentApp
        mutable autoPtr<eddyCurrentApp::Manager> eddyCurrentAppManager_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);


public:

    // Runtime type information
    TypeName("buoyantBoussinesqPimpleEddyCurrentApp::Manager");

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

            //- buoyantBoussinesqPimpleApp
            buoyantBoussinesqPimpleApp::Manager& buoyantBoussinesqPimpleAppManager() const;

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

} // End namespace buoyantBoussinesqPimpleEddyCurrentApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "buoyantBoussinesqPimpleEddyCurrentAppManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
