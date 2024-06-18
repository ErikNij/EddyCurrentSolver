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
    Foam::interApp::Manager

Description
    ...

SourceFiles
    interAppManager.C
    interAppManagerRegion_DEFAULT.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef interAppManager_H
#define interAppManager_H

#include "interAppInclude.H"
#include "pimpleAppControl.H"

#include "solverManagers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interApp
{

/*---------------------------------------------------------------------------*\
                   Class interApp::Manager Declaration
\*---------------------------------------------------------------------------*/

class Manager
:
    public fvSolverManager
{
public:

    // Public typedefs and forward declarations

        typedef pimpleApp::Control Control;
        class Settings;
        class Storage;

    // Solver Manager Construction

        SM_GLOBALCONTROL_ADD_P1(Control, mesh());


        SM_GLOBALSETTINGS_NEW();

        SM_GLOBALSETTINGS_ADD();


        SM_GLOBALSTORAGE_NEW();
        SM_GLOBALSTORAGE_ADD();


        SM_REGION_NEW(DEFAULT);

            SM_REGIONCONTROL_ADD_GLOBALCONTROL();

            SM_REGIONSETTINGS_NEW();

                //- Volume force switch
                SM_SETTINGS_MEMBER_ADD(Switch, volumeForce);

                //- Electrical conductivity switch
                SM_SETTINGS_MEMBER_ADD(Switch, electricalConuctivity);

                //- Use fixed flux-pressure to correct pressure p at
                //  fixedGradient patches
                SM_SETTINGS_MEMBER_ADD(Switch, snGradpFromFlux);

            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                //- Gravitational acceleration
                SM_REGIONITEM_ADD(uniformDimensionedVectorField, g);

                //- Pressure
                SM_REGIONITEM_ADD(volScalarField, p);

                //- Velocity
                SM_REGIONITEM_ADD(volVectorField, U);

                //- Face flux
                SM_REGIONITEM_ADD(surfaceScalarField, phi);

                //- Volume force
                SM_REGIONITEM_ADD(volVectorField, F);

                //- Density
                SM_REGIONITEM_ADD(volScalarField, rho);

                //- Mass flux
                SM_REGIONITEM_ADD(surfaceScalarField, rhoPhi);

                //- Fluid indicator
                SM_REGIONITEM_ADD(volScalarField, alpha1);

                //- Total volume fraction of phase1
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, alpha1Vfrac);

                //- Transport model
                SM_REGIONITEM_ADD(twoPhaseMixture, transport);

                //- Turbulence model
                SM_REGIONITEM_ADD(incompressible::turbulenceModel, turbulence);

                //- Interface model
                SM_REGIONITEM_ADD(interfaceProperties, interface);

                //- Electrical conductivity
                SM_REGIONITEM_ADD(volScalarField, sigma);

            SM_REGIONSTORAGE_ADD();

        SM_REGION_END();


        SM_REGIONS_NEW();

            SM_REGIONS_REGION_ADD(DEFAULT, mesh());

        SM_REGIONS_ADD();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Manager(const Manager&);

        //- Disallow default bitwise assignment
        void operator=(const Manager&);


public:

    // Runtime type information
    TypeName("interApp::Manager");

    // Constructor

        Manager
        (
            const argList& args,
            Time& time,
            fvMesh& mesh,
            bool master = true
        );

    // Destructor

        virtual ~Manager();

    // Member Functions

        // Evolution

            //- Set Courant Number
            virtual bool setCoNum(scalar& CoNum) const;

            //- Next sequence
            virtual void next() const;

            //- Write sequence
            virtual void write() const;

            //- Finalize sequence
            virtual void finalize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interAppManager.C"
#   include "interAppManagerRegion_DEFAULT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
