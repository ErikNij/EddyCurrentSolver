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
    eddyCurrentFoam

Description
    ...

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "eddyCurrentApp.H"

// TODO: Fix regions for settings in solverProeprties!

// TODO: Stop loop if conductivity not given

// TODO: Stop loop if no current present

// TODO: Only do Biot-Savart if mesh is moving or has been changed recently

// TODO: Derived gradient boundary condition for VRe/VIm in conductor region

// TODO: Remove interfaceLabel? Use alias of calculatedFvPatchField or special
//       patchField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

    // Add time selector options for frequencies
    argList::validOptions.insert("freq", "scalarList");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createRegionMeshUninitialized.H"

    using namespace eddyCurrentApp;
    using namespace eddyCurrentApp::Region;

    Manager manager(args, runTime, regionMesh);

    SM_GLOBALREGIONSCOPE(DEFAULT);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalarList freqList;
    globalSettings.freqOpt =
        args.optionReadIfPresent<scalarList>("freq", freqList);

    manager.read();
    manager.init();

    uniformDimensionedScalarField& f0 = storage.f0();
    uniformDimensionedScalarField& omega0 = storage.omega0();

    // Store single value of constant f0 in freqList if no option is given
    if (!globalSettings.freqOpt)
    {
        freqList.setSize(1);
        freqList[0] = f0.value();
    }

    while(manager.selected(freqList))
    {
        // Use time value as frequency value
        f0 = dimensionedScalar(word(), f0.dimensions(), runTime.time().value());
        omega0 = mathematicalConstant::twoPi * f0;

#       include "materialProperties.H"

#       include "A0BiotSavart.H"

#       include "AVInit.H"

#       include "AVLoop.H"

#       include "BUpdate.H"

#       include "derivedFields.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    return(0);
}

// ************************************************************************* //
