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
    interEddyCurrentFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "interEddyCurrentApp.H"

// TODO: Fix emUpdate!

// TODO: Make biot-savart only once but properly!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    interEddyCurrentApp::Manager masterManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        masterManager.eddyCurrentAppManager();

    interApp::Manager& interAppManager =
        masterManager.interAppManager();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    masterManager.read();
    masterManager.init();

    // Init interApp
    {
        using namespace interApp;
        using namespace interApp::Region;

        Manager& manager = interAppManager;

        SM_MANAGERSCOPE();
        SM_REGIONSCOPE(DEFAULT);

#       include "calcSigma.H"
    }

    // Init interEddyCurrentApp (1)
    {
        using namespace interEddyCurrentApp;
        using namespace interEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        eddyCurrentApp::Manager::Storage& ecs =
            eddyCurrentAppManager.storage();

        ecs.sigma().rmap(Region::CONDUCTOR);
        ecs.sigma().rmap(Region::FLUID);
        ecs.sigma().mapInternalField(Region::CONDUCTOR);
        ecs.sigma()[Region::CONDUCTOR].correctBoundaryConditions();
    }

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

    // Init interEddyCurrentApp (2)
    {
        using namespace interEddyCurrentApp;
        using namespace interEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        eddyCurrentApp::Manager::Storage& ecs =
            eddyCurrentAppManager.storage();

        ecs.F().rmap(Region::CONDUCTOR);
        ecs.F().mapInternalField(Region::FLUID);
        fvc::extrapolate(ecs.F()[Region::FLUID]);
    }

// TODO: TEST
    const volScalarField& alpha1 =
        interAppManager.regions().region_DEFAULT().storage().alpha1();
    volScalarField alpha1_old(alpha1);
    Switch emUpdateFirst(true);

    dictionary& pisoDict = interAppManager.regions().region_DEFAULT().mesh().solutionDict().subDict("PISO");
    Switch emUpdateAlpha = pisoDict.lookupOrDefault("emUpdateAlpha", false);
    scalar emUpdateAlphaChange = pisoDict.lookupOrDefault("emUpdateAlphaChange", 0.5);

    while (masterManager.run())
    {
// TODO: TEST
        // Check for magnetic update
        Switch emUpdate(false);
        if (emUpdateAlpha)
        {
            if ((max(mag(alpha1_old - alpha1)).value() > emUpdateAlphaChange) || emUpdateFirst)
            {
                emUpdateFirst = false;
                emUpdate = true;

                alpha1_old = alpha1;
            }
        }
        else
        {
            emUpdate = true;
        }
//         Switch emUpdate(true);
//         Switch emUpdate =
//             eddyCurrentAppManager.control().needsUpdate
//             (
//                 interEddyCurrentApp::Region::FLUID
//             );

        if (emUpdate)
        {
            using namespace interEddyCurrentApp;
            using namespace interEddyCurrentApp::Region;

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

        // Map/Extrapolate and update volume force in fluid region
        if (emUpdate)
        {
            using namespace interEddyCurrentApp;
            using namespace interEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            if (Control::debug)
            {
                Info<< Control::typeName << " : "
                    << "Map/Extrapolate Lorentz-force in fluid region."
                    << endl;
            }

            eddyCurrentApp::Manager::Storage& ecs =
                eddyCurrentAppManager.storage();

            ecs.F().rmap(Region::CONDUCTOR);
            ecs.F().mapInternalField(Region::FLUID);
            fvc::extrapolate(ecs.F()[Region::FLUID]);
        }

        // Solve fluid flow
        {
            using namespace interApp;
            using namespace interApp::Region;

            Manager& manager = interAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "alphaUpLoop.H"
        }

        // Update sigma in conductor region
        {
            using namespace interEddyCurrentApp;
            using namespace interEddyCurrentApp::Region;

            Manager& manager = masterManager;

            SM_MANAGERSCOPE();

            eddyCurrentApp::Manager::Storage& ecs =
                eddyCurrentAppManager.storage();

            ecs.sigma().rmap(Region::CONDUCTOR);
            ecs.sigma().rmap(Region::FLUID);
            ecs.sigma().mapInternalField(Region::CONDUCTOR);
            ecs.sigma()[Region::CONDUCTOR].correctBoundaryConditions();
        }
    }

    return(0);
}


// ************************************************************************* //
