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
    makeWedgeMesh

Description
    Make wedge mesh from a 2D mesh

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validOptions.insert("angle", "scalar [deg]");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Store initial instance of mesh
    word pointsInstance = mesh.pointsInstance();

    // Read optional angle
    scalar angle = 5.0;
    args.optionReadIfPresent("angle", angle);

    scalar phi = mathematicalConstant::pi/180.0 * angle/2.0;

    // Analyse mesh geometry and extract
    label nGeometricD = mesh.nGeometricD();
    Vector<label> geometricD = mesh.geometricD();

    labelHashSet nSet;
    labelHashSet tSet;

    // Only proceed if this mesh is 2D
    if (nGeometricD == 2)
    {
        // Extract normal and tangential directions
        for (direction cmpt = 0; cmpt < Vector<label>::nComponents; cmpt++)
        {
            if (geometricD[cmpt] < 0)
            {
                nSet.insert(cmpt);
            }
            else
            {
                tSet.insert(cmpt);
            }
        }

        labelList n = nSet.sortedToc();
        labelList t = tSet.sortedToc();

        // Create the new points
        pointIOField newPoints
        (
            IOobject
            (
                "points",
                pointsInstance,
                mesh.meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh.allPoints()
        );

        Info << "Transform points" << endl;

        bool mesh2D = false;

        const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

        // Search for empty patches
        forAll (boundaryMesh, patchI)
        {
            const polyPatch& patch = boundaryMesh[patchI];

            if (isA<emptyPolyPatch>(patch))
            {
                const labelList& meshPoints = patch.meshPoints();

                // Transform mesh points
                forAll (meshPoints, patchPointI)
                {
                    label pointI = meshPoints[patchPointI];

                    point p = newPoints[pointI];

                    scalar s = Foam::sign(p[n[0]]);
                    scalar r = p[t[0]];

                    p[n[0]] = s*r*Foam::sin(phi);
                    p[t[0]] =   r*Foam::cos(phi);

                    newPoints[pointI] = p;
                }

                mesh2D = true;
            }
        }

        if (mesh2D)
        {
            Info << "Write new mesh" << endl;

            // Write new points
            newPoints.write();

            // Update mesh according to new points
            mesh.readUpdate();

            // Turn empty patches into wedge ones. This is currently
            // a cruel hack in order to avoid juggling with the all
            // details of boundaryMesh and/or the boundary dictionary itself
            forAll (boundaryMesh, patchI)
            {
                const polyPatch& patch = boundaryMesh[patchI];

                if (isA<emptyPolyPatch>(patch))
                {
                    const_cast<word&>(patch.type()) = "wedge";
                }
            }

            // WARNING: Now the mash became actually invalid, as we have just
            // brutaly RENAMED the empty patches instead of replacing them
            // with real wedgePatches. Due to our laziness, any further usage
            // (after writing the new boundaryMesh) is stronlgy discouraged!

            // Write new patch types
            boundaryMesh.write();

            Info << "Use 'collapseEdges' to remove axis patch!" << endl;
        }
        else
        {
            Info << "Mesh is 2D without empty patches! Already wedge?" << endl;
        }
    }
    else
    {
        Info << "Mesh is not 2D!" << endl;
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
