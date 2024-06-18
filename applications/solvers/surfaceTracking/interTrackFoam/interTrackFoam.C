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
    interTrackFoam

Description

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "interTrackApp.H"

// TODO: Port ddtCorr from OpenFOAM 3.0.x or use new approach from Hrv

// TODO: Fix meshPhi issues!

// TODO: Pressure naming p vs. p_rgh (gh, ghf, hRef, ...)

// TODO: Pressure boundary conditions (F, fixedFluxPressure)

// TODO: Two fluids

// TODO: Two fluids: Pressure equation reference

// TODO: Continuity errors

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    using namespace Foam;

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    using namespace interTrackApp;
    using namespace interTrackApp::Region;

    Manager manager(args, runTime, mesh);

    SM_GLOBALREGIONSCOPE(DEFAULT);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Two fluids: Pressure equation reference
// #   include "setRefCell.H"

// TODO: Continuity errors
// #   include "initContinuityErrs.H"

    while (manager.run())
    {
#       include "meshUpdate.H"

#       include "UpLoop.H"

// TODO: Continuity errors
// #       include "volContinuity.H"
    }

    return(0);
}

// ************************************************************************* //
