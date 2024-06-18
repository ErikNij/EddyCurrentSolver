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
    Foam::buoyantBoussinesqPimpleApp::Manager

Description
    ...

SourceFiles
    buoyantBoussinesqPimpleAppManager.C
    buoyantBoussinesqPimpleAppManagerRegion_DEFAULT.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef buoyantBoussinesqPimpleAppManager_H
#define buoyantBoussinesqPimpleAppManager_H

#include "buoyantBoussinesqPimpleAppInclude.H"
#include "pimpleAppControl.H"

#include "solverManagers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace buoyantBoussinesqPimpleApp
{

/*---------------------------------------------------------------------------*\
             Class buoyantBoussinesqPimpleApp::Manager Declaration
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

                //- Heat source switch
                SM_SETTINGS_MEMBER_ADD(Switch, heatSource);

                //- Use fixed flux-pressure to correct pressure p at
                //  fixedGradient patches
                SM_SETTINGS_MEMBER_ADD(Switch, snGradpFromFlux);

            SM_REGIONSETTINGS_ADD();

            SM_REGIONSTORAGE_NEW();

                //- Gravitational acceleration
                SM_REGIONITEM_ADD(uniformDimensionedVectorField, g);

                //- Temperatur
                SM_REGIONITEM_ADD(volScalarField, T);

                //- Pressure
                SM_REGIONITEM_ADD(volScalarField, p);

                //- Velocity
                SM_REGIONITEM_ADD(volVectorField, U);

                //- Face flux
                SM_REGIONITEM_ADD(surfaceScalarField, phi);

                //- Transport model
                SM_REGIONITEM_ADD(singlePhaseTransportModel, transport);

                //- Thermal expansion coefficient
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, beta);

                //- Reference temperature
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, TRef);

                //- Laminar Prandtl number
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, Pr);

                //- Turbulent Prandtl number
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, Prt);

                //- Turbulence model
                SM_REGIONITEM_ADD(incompressible::turbulenceModel, turbulence);

                //- Kinematic turbulent thermal conductivity
                SM_REGIONITEM_ADD(volScalarField, alphat);

                //- Kinematic density
                SM_REGIONITEM_ADD(volScalarField, rhok);

                //- Volume force
                SM_REGIONITEM_ADD(volVectorField, F);

                //- Reference density
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, rhoRef);

                //- Heat source
                SM_REGIONITEM_ADD(volScalarField, Q);

                //- Reference specific heat capacity
                SM_REGIONITEM_ADD(uniformDimensionedScalarField, CpRef);

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
    TypeName("buoyantBoussinesqPimpleApp::Manager");

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

} // End namespace buoyantBoussinesqPimpleApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "buoyantBoussinesqPimpleAppManager.C"
#   include "buoyantBoussinesqPimpleAppManagerRegion_DEFAULT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
