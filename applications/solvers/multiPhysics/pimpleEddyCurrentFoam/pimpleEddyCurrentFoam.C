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
    pimpleEddyCurrentFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "pimpleEddyCurrentApp.H"

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

    pimpleEddyCurrentApp::Manager masterManager(args, runTime, regionMesh);

    eddyCurrentApp::Manager& eddyCurrentAppManager =
        masterManager.eddyCurrentAppManager();

    pimpleApp::Manager& pimpleAppManager =
        masterManager.pimpleAppManager();

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

    // Init pimpleEddyCurrentApp
    {
        using namespace pimpleEddyCurrentApp;
        using namespace pimpleEddyCurrentApp::Region;

        Manager& manager = masterManager;

        SM_MANAGERSCOPE();

        eddyCurrentApp::Manager::Storage& ecs =
            eddyCurrentAppManager.storage();

        ecs.F().rmap(Region::CONDUCTOR);
        ecs.F().mapInternalField(Region::FLUID);
        fvc::extrapolate(ecs.F()[Region::FLUID]);
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
            using namespace pimpleEddyCurrentApp;
            using namespace pimpleEddyCurrentApp::Region;

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
            using namespace pimpleEddyCurrentApp;
            using namespace pimpleEddyCurrentApp::Region;

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
            using namespace pimpleApp;
            using namespace pimpleApp::Region;

            Manager& manager = pimpleAppManager;

            SM_MANAGERSCOPE();
            SM_REGIONSCOPE(DEFAULT);

#           include "UpLoop.H"
        }
    }

    return(0);
}


// ************************************************************************* //
