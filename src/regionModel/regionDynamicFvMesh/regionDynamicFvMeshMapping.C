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

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "regionDynamicFvMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionDynamicFvMesh::patchMapMeshVelocityDirectMapped
(
    label fromRegionI,
    label toRegionI
) const
{
    if (isFeMotionSolver(toRegionI))
    {
        tetPointVectorField& toMotionU =
            const_cast<tetPointVectorField&>
            (
                operator[](toRegionI).objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        labelList patchMap =
            regionPolyMesh::patchMapDirectMapped
            (
                fromRegionI,
                toRegionI
            )[0];

        forAll (operator[](fromRegionI).boundary(), fromPatchI)
        {
            label toPatchI = patchMap[fromPatchI];

            if (toPatchI != -1)
            {
                const polyPatch& fromPolyPatch =
                    operator[](fromRegionI).boundaryMesh()[fromPatchI];

                const polyPatch& toPolyPatch =
                    operator[](toRegionI).boundaryMesh()[toPatchI];

                tmp<vectorField> ttoTotalDisplacement
                (
                    new vectorField
                    (
                        fromPolyPatch.localPoints()
                      - toPolyPatch.localPoints()
                    )
                );

                vectorField& toTotalDisplacement = ttoTotalDisplacement();

                fixedValueTetPolyPatchVectorField& toMotionUpatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        toMotionU.boundaryField()[toPatchI]
                    );

                tetPolyPatchInterpolation toTppiPatch
                (
                    refCast<const faceTetPolyPatch>
                    (
                        toMotionUpatch.patch()
                    )
                );

                toMotionUpatch ==
                    toTppiPatch.pointToPointInterpolate
                    (
                        toTotalDisplacement/time().deltaT().value()
                    );

                ttoTotalDisplacement.clear();
            }
        }

        toMotionU.correctBoundaryConditions();
    }
    else if (isFvMotionSolver(toRegionI))
    {
        // Mesh velocity field
        pointVectorField& toPointMotionU =
            const_cast<pointVectorField&>
            (
                operator[](toRegionI).objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        labelList patchMap =
            patchMapDirectMapped
            (
                fromRegionI,
                toRegionI
            )[0];

        forAll (operator[](fromRegionI).boundary(), fromPatchI)
        {
            label toPatchI = patchMap[fromPatchI];

            if (toPatchI != -1)
            {
                const polyPatch& fromPolyPatch =
                    operator[](fromRegionI).boundaryMesh()[fromPatchI];

                const polyPatch& toPolyPatch =
                    operator[](toRegionI).boundaryMesh()[toPatchI];

                // WARNING: Assuming same point/face/cell ordering on both
                //          sides of this direct mapped patch!

                tmp<vectorField> ttoTotalDisplacement
                (
                    new vectorField
                    (
                        fromPolyPatch.localPoints()
                      - toPolyPatch.localPoints()
                    )
                );

                vectorField& toTotalDisplacement = ttoTotalDisplacement();

                fixedValuePointPatchVectorField& toPointMotionUpatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        toPointMotionU.boundaryField()[toPatchI]
                    );

                toPointMotionUpatch ==
                    toTotalDisplacement/time().deltaT().value();

                ttoTotalDisplacement.clear();
            }
        }

        toPointMotionU.correctBoundaryConditions();
    }
}


void Foam::regionDynamicFvMesh::patchMapMeshVelocityDirectMapped
(
    const word& fromRegionName,
    const word& toRegionName
) const
{
    label fromRegionI = regions()[fromRegionName];
    label toRegionI = regions()[toRegionName];

    return patchMapMeshVelocityDirectMapped
    (
        fromRegionI,
        toRegionI
    );
}


void Foam::regionDynamicFvMesh::patchMapMeshPointsDirectMapped
(
    label fromRegionI,
    label toRegionI,
    pointField& toPoints
) const
{
    if (toPoints.size() == 0)
    {
        toPoints = operator[](toRegionI).allPoints();
    }

    labelList patchMap =
        patchMapDirectMapped
        (
            fromRegionI,
            toRegionI
        )[0];

    forAll (operator[](fromRegionI).boundary(), fromPatchI)
    {
        label toPatchI = patchMap[fromPatchI];

        if (toPatchI != -1)
        {
            const polyPatch& fromPolyPatch =
                operator[](fromRegionI).boundaryMesh()[fromPatchI];

            const polyPatch& toPolyPatch =
                operator[](toRegionI).boundaryMesh()[toPatchI];

            // WARNING: Assuming same point/face/cell ordering on both
            //          sides of this direct mapped patch!

            forAll (fromPolyPatch.meshPoints(), fromPatchPointI)
            {
                label toPointI = toPolyPatch.meshPoints()[fromPatchPointI];

                toPoints[toPointI] = fromPolyPatch.localPoints()[fromPatchPointI];
            }
        }
    }
}


void Foam::regionDynamicFvMesh::patchMapMeshPointsDirectMapped
(
    const word& fromRegionName,
    const word& toRegionName,
    pointField& toPoints
) const
{
    label fromRegionI = regions()[fromRegionName];
    label toRegionI = regions()[toRegionName];

    return patchMapMeshPointsDirectMapped
    (
        fromRegionI,
        toRegionI,
        toPoints
    );
}


// ************************************************************************* //

