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

\*---------------------------------------------------------------------------*/

#include "faSubMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faSubMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<pointField> faSubMesh::calcNewPoints() const
{
    const int nSubPolyMeshPoints = subPolyMesh().points().size();
    const int nSubAreaMeshPoints = subAreaMesh().points().size();

    const int nUsedBasePolyMeshPoints =
        nSubPolyMeshPoints - nSubAreaMeshPoints;
    const int nUsedBaseAreaMeshPoints = baseAreaMesh().points().size();

    const pointField& basePolyMeshPoints = basePolyMesh().points();
    const pointField& baseAreaMeshPoints = baseAreaMesh().points();
    const pointField& splitPoints = subSplitPoints();

    tmp<pointField> tpoints
    (
        new pointField(subPolyMesh().points())
    );
    pointField& points = tpoints();

    // Used vertices from base polyMesh come first
    for (label pointI = 0; pointI < nUsedBasePolyMeshPoints; pointI++)
    {
        points[pointI] =
            basePolyMeshPoints[pointSubToBaseMap()[pointI]];
    }

    // Used vertices from base area mesh come second
    for (label pointI = 0; pointI < nUsedBaseAreaMeshPoints; pointI++)
    {
        // Point map is not necessary here, as points are ordered identically
        points[nUsedBasePolyMeshPoints + pointI] =
            baseAreaMeshPoints[pointI];
    }

    // Split points come last
    forAll(splitPoints, pointI)
    {
        // Point map is not necessary here, as points are ordered identically
        points[nUsedBasePolyMeshPoints + nUsedBaseAreaMeshPoints + pointI] =
            splitPoints[pointI];
    }

    return tpoints;
}


void faSubMesh::clearGeom() const
{
    deleteDemandDrivenData(faceCurvaturesPtr_);
}


void faSubMesh::clearOut() const
{
    clearGeom();

    deleteDemandDrivenData(faceSubToBaseAreaMapPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faSubMesh::faSubMesh
(
    const faMesh& baseAreaMesh,
    const pointField& subSplitPoints
)
:
basePolyMesh_(baseAreaMesh.mesh()),
baseAreaMesh_(baseAreaMesh),
subPolyMesh_
(
    IOobject
    (
        basePolyMesh_.name() + "_faSubMesh",
        basePolyMesh_.time().timeName(),
        basePolyMesh_.time(),
        IOobject::MUST_READ
    )
),
subToBaseAddressing_(subPolyMesh_),
subAreaMesh_(subPolyMesh_),
subSplitPoints_(subSplitPoints),
faceSubToBaseAreaMapPtr_(NULL),
faceCurvaturesPtr_(NULL)
{
    subAreaMesh_.correctPatchPointNormals() =
        baseAreaMesh_.correctPatchPointNormals();
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

faSubMesh::~faSubMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
