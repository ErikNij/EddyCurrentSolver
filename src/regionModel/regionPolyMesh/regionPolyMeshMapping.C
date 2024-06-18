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

#include "directMappedPatchBase.H"
#include "directMappedPolyPatch.H"
#include "directMappedWallPolyPatch.H"

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regionPolyMesh::map
(
    pointField& givenPoints,
    label regionI,
    const word& patchName
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    const pointField& points0 = operator[](regionI0).points();
    const pointField& points = operator[](regionI).points();

    const labelList& map =
        typeMap(addressingTypes::POINT, regionI, regionI0);

    if (patchName == "")
    {
        forAll (points, pointI)
        {
            givenPoints[pointI] = points0[map[pointI]];
        }
    }
    else
    {
        const polyBoundaryMesh& boundaryMesh =
            operator[](regionI).boundaryMesh();

        label patchI = boundaryMesh.findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn("regionPolyMesh::rmap(...)")
                << "Given patch name " << patchName
                << " for point mapping does not exist in"
                << " region " << regions()[regionI]
                << abort(FatalError);
        }

        const polyPatch& patch = boundaryMesh[patchI];

        const labelList& meshPoints = patch.meshPoints();

        forAll (meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];

            givenPoints[pointI] = points[map[pointI]];
        }
    }
}


void regionPolyMesh::map
(
    pointField& givenPoints,
    const word& regionName,
    const word& patchName
) const
{
    map
    (
        givenPoints,
        regions()[regionName],
        patchName
    );
}


tmp<pointField> regionPolyMesh::map
(
    label regionI,
    const word& patchName
) const
{
    tmp<pointField> tNewPoints
    (
        new pointField(operator[](regionI).allPoints())
    );

    pointField& newPoints = tNewPoints();

    map
    (
        newPoints,
        regionI,
        patchName
    );

    return tNewPoints;
}


tmp<pointField> regionPolyMesh::map
(
    const word& regionName,
    const word& patchName
) const
{
    return map
    (
        regions()[regionName],
        patchName
    );
}


void regionPolyMesh::rmap
(
    pointField& givenPoints,
    label regionI,
    const word& patchName
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    const pointField& points = operator[](regionI).points();

    const labelList& map =
        typeMap(addressingTypes::POINT, regionI, regionI0);

    if (patchName == "")
    {
        forAll (points, pointI)
        {
            givenPoints[map[pointI]] = points[pointI];
        }
    }
    else
    {
        const polyBoundaryMesh& boundaryMesh =
            operator[](regionI).boundaryMesh();

        label patchI = boundaryMesh.findPatchID(patchName);

        if (patchI == -1)
        {
            FatalErrorIn("regionPolyMesh::rmap(...)")
                << "Given patch name " << patchName
                << " for point mapping does not exist in"
                << " region " << regions()[regionI]
                << abort(FatalError);
        }

        const polyPatch& patch = boundaryMesh[patchI];

        const labelList& meshPoints = patch.meshPoints();

        forAll (meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];

            givenPoints[map[pointI]] = points[pointI];
        }
    }
}


void regionPolyMesh::rmap
(
    pointField& givenPoints,
    const word& regionName,
    const word& patchName
) const
{
    rmap
    (
        givenPoints,
        regions()[regionName],
        patchName
    );
}


tmp<pointField> regionPolyMesh::rmap
(
    label regionI,
    const word& patchName
) const
{
    label regionI0 = regions()[polyMesh::defaultRegion];

    tmp<pointField> tNewPoints
    (
        new pointField(operator[](regionI0).allPoints())
    );

    pointField& newPoints = tNewPoints();

    rmap
    (
        newPoints,
        regionI,
        patchName
    );

    return tNewPoints;
}


tmp<pointField> regionPolyMesh::rmap
(
    const word& regionName,
    const word& patchName
) const
{
    return rmap
    (
        regions()[regionName],
        patchName
    );
}


labelListList regionPolyMesh::patchMapDirectMapped
(
    label fromRegionI,
    label toRegionI
) const
{
    const polyMesh& fromMesh = operator[](fromRegionI);
    const polyMesh& toMesh = operator[](toRegionI);

    label pmapSize = fromMesh.boundaryMesh().size();

    labelListList pmap
    (
        2,
        labelList
        (
            pmapSize,
            -1
        )
    );

    forAll (fromMesh.boundaryMesh(), fromPatchI)
    {
        if
        (
            (
                fromMesh.boundaryMesh()[fromPatchI].type()
            == directMappedPolyPatch::typeName
            )
            ||
            (
                fromMesh.boundaryMesh()[fromPatchI].type()
            == directMappedWallPolyPatch::typeName
            )
        )
        {
            const directMappedPatchBase& dmpb =
                refCast<const directMappedPatchBase>
                (
                    fromMesh.boundaryMesh()[fromPatchI]
                );

            word toPatchName = dmpb.samplePatch();

            label toPatchI =
                toMesh.boundaryMesh().findPatchID(toPatchName);

            pmap[0][fromPatchI] = toPatchI;
            pmap[1][toPatchI] = fromPatchI;
        }
    }

    return pmap;
}


labelListList regionPolyMesh::patchMapDirectMapped
(
    const word& fromRegionName,
    const word& toRegionName
) const
{
    label fromRegionI = regions()[fromRegionName];
    label toRegionI = regions()[toRegionName];

    return patchMapDirectMapped
    (
        fromRegionI,
        toRegionI
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

