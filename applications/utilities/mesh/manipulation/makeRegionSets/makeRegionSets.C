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
    makeRegionSets

Description
    Make sets for regions.

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "meshTools.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"

// TODO: Warning/Exit if region option is not used (for defaultRegion). It
//       simply does not make any sense.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

    // Add option to write empty sets
    argList::validOptions.insert("writeEmptySets", "");
    argList::validOptions.insert("liveObjectsOnly", "");

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

#   include "createNamedPolyMesh.H"

    bool writeEmptySets = args.optionFound("writeEmptySets");
    bool liveObjectsOnly = args.optionFound("liveObjectsOnly");

    const word facesInstance =
        runTime.findInstance(polyMesh::meshSubDir, "faces");


    // Search for list of objects for the time of default region
    IOobjectList objects
    (
        runTime,
        facesInstance,
        polyMesh::meshSubDir/"sets"
    );

    Info<< "Searched : " << facesInstance/polyMesh::meshSubDir/"sets"
        << nl
        << "Found    : " << objects.names() << nl
        << endl;


    // pointSets
    IOobjectList pointObjects(objects.lookupClass(pointSet::typeName));
    forAllConstIter(IOobjectList, pointObjects, iter)
    {
        pointSet set(*iter());

        IOobject pointRegionAddObj
        (
            "pointRegionAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        IOobject pointMapObj
        (
            "pointMap",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (!pointRegionAddObj.headerOk() && pointMapObj.headerOk())
        {
            pointRegionAddObj = pointMapObj;
        }

        pointRegionAddObj.readOpt() = IOobject::MUST_READ;

        labelIOList pointRegionAdd = labelIOList(pointRegionAddObj);

        const label nRegionPoints = mesh.nPoints();

        Info<< "Split pointSet " << set.name() << " ... ";

        labelHashSet regionSet;

        forAll (pointRegionAdd, pointI)
        {
            if (liveObjectsOnly && pointI >= nRegionPoints)
            {
                break;
            }

            if (set.found(pointRegionAdd[pointI]))
            {
                regionSet.insert(pointI);
            }
        }

        if (!regionSet.empty() || writeEmptySets)
        {
            pointSet ps
            (
                mesh,
                set.name(),
                regionSet,
                IOobject::NO_WRITE
            );

            ps.write();
        }

        Info<< "OK" << endl;
    }


    // faceSets
    IOobjectList faceObjects(objects.lookupClass(faceSet::typeName));
    forAllConstIter(IOobjectList, faceObjects, iter)
    {
        // Set not in memory. Load it.
        faceSet set(*iter());

        IOobject faceRegionAddObj
        (
            "faceRegionAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        IOobject faceMapObj
        (
            "faceMap",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        bool fromMap = false;

        if (!faceRegionAddObj.headerOk() && faceMapObj.headerOk())
        {
            faceRegionAddObj = faceMapObj;

            fromMap = true;
        }

        faceRegionAddObj.readOpt() = IOobject::MUST_READ;

        labelIOList faceRegionAdd = labelIOList(faceRegionAddObj);

        const label nRegionFaces = mesh.nFaces();

        Info<< "Split faceSet " << set.name() << " ... ";

        labelHashSet regionSet;

        forAll (faceRegionAdd, faceI)
        {
            if (liveObjectsOnly && faceI >= nRegionFaces)
            {
                break;
            }

            if (!fromMap)
            {
                if (set.found(mag(faceRegionAdd[faceI]) - 1))
                {
                    regionSet.insert(faceI);
                }
            }
            else
            {
                if (set.found(faceRegionAdd[faceI]))
                {
                    regionSet.insert(faceI);
                }
            }
        }

        if (!regionSet.empty() || writeEmptySets)
        {
            faceSet fs
            (
                mesh,
                set.name(),
                regionSet,
                IOobject::NO_WRITE
            );

            fs.write();
        }

        Info<< "OK" << endl;
    }


    // cellSets
    IOobjectList cellObjects(objects.lookupClass(cellSet::typeName));
    forAllConstIter(IOobjectList, cellObjects, iter)
    {
        cellSet set(*iter());

        IOobject cellRegionAddObj
        (
            "cellRegionAddressing",
            mesh.facesInstance(),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        IOobject cellMapObj
        (
            "cellMap",
            mesh.facesInstance(),
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

        const label nRegionCells = mesh.nCells();

        Info<< "Split cellSet " << set.name() << " ... ";

        labelHashSet regionSet;

        forAll (cellRegionAdd, cellI)
        {
            if (liveObjectsOnly && cellI >= nRegionCells)
            {
                break;
            }

            if (set.found(cellRegionAdd[cellI]))
            {
                regionSet.insert(cellI);
            }
        }

        if (!regionSet.empty() || writeEmptySets)
        {
            cellSet cs
            (
                mesh,
                set.name(),
                regionSet,
                IOobject::NO_WRITE
            );

            cs.write();
        }

        Info<< "OK" << endl;
    }


    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
