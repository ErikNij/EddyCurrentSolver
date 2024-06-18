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
    makeRegionCellDecomposition

Description
    Make cellDecomposition for regions.

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "meshTools.H"

// TODO: Warning/Exit if region option is not used (for defaultRegion). It
//       simply does not make any sense.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            runTime.findInstance(polyMesh::meshSubDir, "faces"),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

#   include "createNamedPolyMesh.H"

    const word facesInstance = mesh.facesInstance();

    IOobject cellRegionAddObj
    (
        "cellRegionAddressing",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    IOobject cellMapObj
    (
        "cellMap",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!cellRegionAddObj.headerOk() && cellMapObj.headerOk())
    {
        cellRegionAddObj = cellMapObj;
    }

    cellRegionAddObj.readOpt() = IOobject::MUST_READ;

    labelIOList cellRegionAdd = labelIOList(cellRegionAddObj);

    Info<< "Split cellDecomposition ... ";

    labelIOList regionCellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            facesInstance,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh.cells().size())
    );

    forAll (regionCellDecomposition, cellI)
    {
        regionCellDecomposition[cellI] =
            cellDecomposition[cellRegionAdd[cellI]];
    }

    regionCellDecomposition.write();

    Info<< "OK" << endl;

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
