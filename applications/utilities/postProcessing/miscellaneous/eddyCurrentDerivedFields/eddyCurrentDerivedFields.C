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
    eddyCurrentDerivedFields

Description
    Calculate derived fields of eddyCurrentFoam

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "eddyCurrentApp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

    using namespace eddyCurrentApp;
    using namespace eddyCurrentApp::Region;

    Manager manager(args, runTime, regionMesh);

    SM_GLOBALREGIONSCOPE(DEFAULT);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    manager.read();

    globalSettings.lorentzForce = true;
    globalSettings.magneticPressure = true;
    globalSettings.jouleHeat = true;

    manager.init();

    uniformDimensionedScalarField& omega0 = storage.omega0();

    while(manager.once())
    {
#       include "materialProperties.H"

#       include "AVInit.H"

#       include "BUpdate.H"

#       include "derivedFields.H"

        globalStorage.ARe()[CONDUCTOR].write();
        globalStorage.AIm()[CONDUCTOR].write();

        globalStorage.BRe()[CONDUCTOR].write();
        globalStorage.BIm()[CONDUCTOR].write();

        globalStorage.VReGrad()[CONDUCTOR].write();
        globalStorage.VImGrad()[CONDUCTOR].write();

        manager.regions().region_CONDUCTOR().storage().jRe().write();
        manager.regions().region_CONDUCTOR().storage().jIm().write();

        globalStorage.F()[CONDUCTOR].write();
        globalStorage.pB()[CONDUCTOR].write();
        globalStorage.Q()[CONDUCTOR].write();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    return(0);
}

// ************************************************************************* //
