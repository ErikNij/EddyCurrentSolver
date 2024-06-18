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
    zonesToSets

Description
    Convert zones to sets with topSet-specific prefix.

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "meshTools.H"

#include "cellZone.H"
#include "faceZone.H"
#include "pointZone.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedPolyMesh.H"



    // cellZones

    Info << "cellZones:" << endl
        << mesh.cellZones().size() << endl
        << "(" << endl;

    forAll(mesh.cellZones(), cellZoneI)
    {
        cellZone& curCellZone = mesh.cellZones()[cellZoneI];

        Info << "    " << curCellZone.name() << endl;

        labelHashSet curCellZoneLabelHashSet(curCellZone);

        cellSet curCellSet(mesh, "cellSet_"+curCellZone.name(), curCellZoneLabelHashSet);

        curCellSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    // faceZones

    Info << "faceZones:" << endl
        << mesh.faceZones().size() << endl
        << "(" << endl;

    forAll(mesh.faceZones(), faceZoneI)
    {
        faceZone& curFaceZone = mesh.faceZones()[faceZoneI];

        Info << "    " << curFaceZone.name() << endl;

        labelHashSet curFaceZoneLabelHashSet(curFaceZone);

        faceSet curFaceSet(mesh, "faceSet_"+curFaceZone.name(), curFaceZoneLabelHashSet);

        curFaceSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    // pointZones

    Info << "pointZones:" << endl
        << mesh.pointZones().size() << endl
        << "(" << endl;

    forAll(mesh.pointZones(), pointZoneI)
    {
        pointZone& curPointZone = mesh.pointZones()[pointZoneI];

        Info << "    " << curPointZone.name() << endl;

        labelHashSet curPointZoneLabelHashSet(curPointZone);

        pointSet curPointSet(mesh, "pointSet_"+curPointZone.name(), curPointZoneLabelHashSet);

        curPointSet.Foam::regIOobject::write();
    }

    Info << ")" << endl << nl;



    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
