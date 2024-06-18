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

Application
    buoyantBoussinesqPimpleEddyCurrentFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "buoyantBoussinesqPimpleEddyCurrentApp.H"

// TODO: Fix emUpdate!

// TODO: Make biot-savart only once but properly!

// TODO: Do not write U/phi

// TODO: Rename to bbqPimpleThermalEddyCurrentFoam

// TODO: Pressure naming p vs. p_rgh (gh, ghf, hRef, ...)

// TODO: Relaxation of pressure?

// TODO: Pressure boundary conditions (F, bouyantPressure vs. fixedFluxPressure)

// TODO: Continuity errors

// TODO: Wall functions for alphat?

// TODO: Radiation (radiation.ST, rhoCpRef, ...)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    buoyantBoussinesqPimpleEddyCurrentApp::Manager
        masterManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        masterManager.eddyCurrentAppManager();

    buoyantBoussinesqPimpleApp::Manager& buoyantBoussinesqPimpleAppManager =
        masterManager.buoyantBoussinesqPimpleAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    masterManager.read();
    masterManager.init();

    // Init eddyCurrentApp
    {
        using namespace eddyCurrentApp;
        using namespace eddyCurrentApp::Region;

        Manager& manager = eddyCurrentAppManager;

        SM_GLOBALREGIONSCOPE(DEFAULT);

        uniformDimensionedScalarField& omega0 = storage.omega0();

        {
#           include "materialProperties.H"

#           include "A0BiotSavart.H"

#           include "AVInit.H"

#           include "BUpdate.H"

#           include "derivedFields.H"
        }
    }

    // Init buoyantBoussinesqPimpleEddyCurrentApp
    {
        using namespace buoyantBoussinesqPimpleEddyCurrentApp;
        using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        Manager::Storage& mms = masterManager.storage();

        eddyCurrentApp::Manager::Storage& ecs =
            eddyCurrentAppManager.storage();

        buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage& bpDs =
            buoyantBoussinesqPimpleAppManager.regions().region_DEFAULT().storage();

        ecs.F().rmap(Region::CONDUCTOR);
        ecs.F().mapInternalField(Region::FLUID);
        fvc::extrapolate(ecs.F()[Region::FLUID]);

        ecs.Q().rmap(Region::CONDUCTOR);
        ecs.Q().mapInternalField(Region::THERMAL);
        fvc::extrapolate(ecs.Q()[Region::THERMAL]);

        mms.rhoCp().rmap(Region::THERMAL);
        mms.rhoCp()[Region::FLUID] = bpDs.rhoRef() * bpDs.CpRef();
        mms.rhoCp().rmap(Region::FLUID);
        mms.rhoCp().mapInternalField(Region::THERMAL);
        fvc::extrapolate(mms.rhoCp()[Region::THERMAL]);

        mms.lambda().rmap(Region::THERMAL);
// TODO: alphat (Boundary conditions)
        mms.lambda()[Region::FLUID] = mms.rhoCp()[Region::FLUID]
                                    * ( bpDs.transport().nu()/bpDs.Pr()
                                      + bpDs.alphat() );
        mms.lambda().rmap(Region::FLUID);
        mms.lambda().mapInternalField(Region::THERMAL);
        fvc::extrapolate(mms.lambda()[Region::THERMAL]);
    }

    while (masterManager.run())
    {
        // Check for magnetic update
        Switch emUpdate(true);
//         Switch emUpdate =
//             eddyCurrentAppManager.control().needsUpdate
//             (
//                 interEddyCurrentApp::Region::FLUID
//             );

        if (emUpdate)
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info<< Control::typeName << " : "
                    << "Update of electromagntic fields due."
                    << endl;
            }
        }

        // Solve eddy-current problem
        if (emUpdate)
        {
            using namespace eddyCurrentApp;
            using namespace eddyCurrentApp::Region;

            Manager& manager = eddyCurrentAppManager;

            SM_GLOBALREGIONSCOPE(DEFAULT);

            uniformDimensionedScalarField& omega0 = storage.omega0();

            {
#               include "materialProperties.H"

#               include "AVLoop.H"

#               include "BUpdate.H"

#               include "derivedFields.H"
            }
        }

        // Map/Extrapolate and update volume force in fluid region and
        // Joule heat in conducting region
        if (emUpdate)
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            eddyCurrentApp::Manager::Storage& ecs =
                eddyCurrentAppManager.storage();

            if (Control::debug)
            {
                Info<< Control::typeName << " : "
                    << "Map/Extrapolate Lorentz-force to fluid region."
                    << endl;
            }

            ecs.F().rmap(Region::CONDUCTOR);
            ecs.F().mapInternalField(Region::FLUID);
            fvc::extrapolate(ecs.F()[Region::FLUID]);

            if (Control::debug)
            {
                Info<< Control::typeName << " : "
                    << "Map/Extrapolate Joule-heat to thermal region."
                    << endl;
            }

            ecs.Q().rmap(Region::CONDUCTOR);
            ecs.Q().mapInternalField(Region::THERMAL);
            fvc::extrapolate(ecs.Q()[Region::THERMAL]);
        }

        // Map/Extrapolate lambda
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            Manager::Storage& mms = masterManager.storage();

            buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage& bpDs =
                buoyantBoussinesqPimpleAppManager.regions().region_DEFAULT().storage();

// TODO: This is constant and is only necessary during init
//             mms.rhoCp().rmap(Region::THERMAL);
//             mms.rhoCp()[Region::FLUID] = bpDs.rhoRef() * bpDs.CpRef();
//             mms.rhoCp().rmap(Region::FLUID);
//             mms.rhoCp().mapInternalField(Region::THERMAL);
//             fvc::extrapolate(mms.rhoCp()[Region::THERMAL]);

            mms.lambda().rmap(Region::THERMAL);
// TODO: alphat (Boundary conditions)
            mms.lambda()[Region::FLUID] = mms.rhoCp()[Region::FLUID]
                                        * ( bpDs.transport().nu()/bpDs.Pr()
                                          + bpDs.alphat() );
            mms.lambda().rmap(Region::FLUID);
            mms.lambda().mapInternalField(Region::THERMAL);
            fvc::extrapolate(mms.lambda()[Region::THERMAL]);
        }

// TODO: Loop over U -> alphat -> T -> rhok -> p?
        // Solve thermal problem
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            Manager::Storage& mms = masterManager.storage();

            eddyCurrentApp::Manager::Storage& ecs =
                eddyCurrentAppManager.storage();

            volScalarField& Q = ecs.Q()[Region::THERMAL];

            volScalarField& T = mms.T()[Region::THERMAL];
// TODO: lambda (Boundary conditions)
            volScalarField& lambda = mms.lambda()[Region::THERMAL];
// TODO: rhoCp (Boundary conditions)
            volScalarField& rhoCp = mms.rhoCp()[Region::THERMAL];

// TODO: Substitute U with phi
// TODO: phi (Boundary conditions)
//             surfaceScalarField& phi = mms.phi()[Region::THERMAL];
            mms.U().rmap(Region::FLUID);
            mms.U().mapInternalField(Region::THERMAL);
            surfaceScalarField phi
            (
                "phi",
                fvc::interpolate(mms.U()[Region::THERMAL])
              & globalMesh[Region::THERMAL].Sf()
            );

            fvScalarMatrix TEqn
            (
                rhoCp*fvm::ddt(T)
              + rhoCp*fvm::div(phi, T)
              - fvm::laplacian(lambda, T)
              - Q
//              ==
//                 rhoCp*radiation.ST(rhoCp, T)
            );

            TEqn.relax();

            TEqn.solve();

//             radiation.correct();
        }

        // Map/Extrapolate and update temperature in fluid region
        {
            using namespace buoyantBoussinesqPimpleEddyCurrentApp;
            using namespace buoyantBoussinesqPimpleEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            Manager::Storage& mms = masterManager.storage();

            if (Control::debug)
            {
                Info<< Control::typeName << " : "
                    << "Map/Extrapolate temperature to fluid region."
                    << endl;
            }

            mms.T().rmap(Region::THERMAL);
            mms.T().mapInternalField(Region::FLUID);
            fvc::extrapolate(mms.T()[Region::FLUID]);
        }

        // Solve fluid flow
        {
            using namespace buoyantBoussinesqPimpleApp;
            using namespace buoyantBoussinesqPimpleApp::Region;

            Manager& manager = buoyantBoussinesqPimpleAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

            if (Control::debug)
            {
                Info<< Control::typeName << " | UpLoop.H : "
                    << "Commencing PIMPLE U-p loop."
                    << endl;
            }

            // --- PIMPLE corrector loop
            while (control.loop())
            {
                uniformDimensionedVectorField& g = storage.g();
                volScalarField& T = storage.T();
                volScalarField& p = storage.p();
                volVectorField& U = storage.U();
                surfaceScalarField& phi = storage.phi();
                uniformDimensionedScalarField& beta = storage.beta();
                uniformDimensionedScalarField& TRef = storage.TRef();
                incompressible::turbulenceModel& turbulence = storage.turbulence();
                volScalarField& rhok = storage.rhok();

                uniformDimensionedScalarField& Prt = storage.Prt();
                volScalarField& alphat = storage.alphat();

#               include "UTpLoop_UEqn.H"
// TODO: alphat wall functions?
// TODO: alphat (Boundary conditions)
#               include "UTpLoop_alphatUpdate.H"
// TODO: Loop over U -> alphat -> T -> rhok -> p?
#               include "UTpLoop_rhokUpdate.H"

                // --- Pressure corrector loop
                while (control.correct())
                {
#                   include "UTpLoop_pEqn.H"
                }

                if (control.turbCorr())
                {
                    turbulence.correct();
                }
            }
        }
    }

    return(0);
}


// ************************************************************************* //
