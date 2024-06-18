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
    Make sets from all cells not in any set.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

    // Add optional argument to specify cellset name
    argList::validOptions.insert("name", "word");

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

#   include "createNamedPolyMesh.H"

    // Read cellSet name
    Foam::word csName = "complementaryCells";
    args.optionReadIfPresent("name", csName);

    // Search for list of objects for the time of the mesh
    IOobjectList objects
    (
        mesh,
        mesh.facesInstance(),
        polyMesh::meshSubDir/"sets"
    );

    Info<< "Searched : " << mesh.facesInstance()/polyMesh::meshSubDir/"sets"
        << nl
        << "Found    : " << objects.names() << nl
        << endl;

    labelHashSet cellsInAnySet;

    IOobjectList cellObjects(objects.lookupClass(cellSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = cellObjects.begin();
        iter != cellObjects.end();
        ++iter
    )
    {
        cellSet set(*iter());

        forAllConstIter (cellSet, set, iter)
        {
            cellsInAnySet.insert(iter.key());
        }
    }

    labelHashSet cellsNotInAnySet;

    forAll (mesh.cells(), cellI)
    {
        if (!cellsInAnySet.found(cellI))
        {
            cellsNotInAnySet.insert(cellI);
        }
    }

    if (!cellsNotInAnySet.empty())
    {
        Info<<"Found " << cellsNotInAnySet.toc().size()
            << " complementary cells. Writing them to cellSet"
            << " '" << csName << "'." << endl;

        cellSet cs
        (
            mesh,
            csName,
            cellsNotInAnySet,
            IOobject::NO_WRITE
        );

        cs.write();
    }
    else
    {
        Info<<"Did not find any complementary cells." << endl;
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //

