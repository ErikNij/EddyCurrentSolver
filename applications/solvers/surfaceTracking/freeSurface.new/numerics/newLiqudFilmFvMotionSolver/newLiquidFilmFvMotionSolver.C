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

#include "newLiquidFilmFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(newLiquidFilmFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        newLiquidFilmFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newLiquidFilmFvMotionSolver::newLiquidFilmFvMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    fvMotionSolver(mesh),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_)
    ),
    pointPoint_(fvMesh_.nPoints(), -1),
    pointPoint2_(fvMesh_.nPoints(), -1),
    max_(0),
    min_(0),
    motionDirection_(lookup("motionDirection"))
{
    motionDirection_ /= mag(motionDirection_) + SMALL;

    boundBox box(mesh.points());

    max_ = (motionDirection_&box.max());
    min_ = (motionDirection_&box.min());

    // Detect the free surface patch

    label aPatchID = fvMesh_.boundaryMesh().findPatchID("freeSurface");

    if(aPatchID == -1)
    {
        FatalErrorIn("newLiquidFilmFvMotionSolver()")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }

    label bottomPatchID =
        fvMesh_.boundaryMesh().findPatchID("bottom");

    if(bottomPatchID == -1)
    {
        FatalErrorIn("newLiquidFilmFvMotionSolver")
            << "Bottom patch not defined.  Please make sure that "
                << " the bottom patches is named as bottom"
                << abort(FatalError);
    }

    const vectorField& points = fvMesh_.points();

    const labelList& meshPointsA =
        fvMesh_.boundaryMesh()[aPatchID].meshPoints();

    forAll(pointPoint_, pointI)
    {
        scalar minDist = GREAT;

        forAll(meshPointsA, mpI)
        {
            vector R = points[meshPointsA[mpI]] - points[pointI];
            R -= motionDirection_*(motionDirection_&R);
            scalar curDist = mag(R);

            if (curDist < minDist)
            {
                minDist = curDist;
                pointPoint_[pointI] = meshPointsA[mpI];
            }
        }
    }


    const labelList& meshPointsBottom =
        fvMesh_.boundaryMesh()[bottomPatchID].meshPoints();

    forAll(pointPoint2_, pointI)
    {
        scalar minDist = GREAT;

        forAll(meshPointsBottom, mpI)
        {
            vector R = points[meshPointsBottom[mpI]] - points[pointI];
            R -= motionDirection_*(motionDirection_&R);
            scalar curDist = mag(R);

            if (curDist < minDist)
            {
                minDist = curDist;
                pointPoint2_[pointI] = meshPointsBottom[mpI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::newLiquidFilmFvMotionSolver::~newLiquidFilmFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::newLiquidFilmFvMotionSolver::curPoints() const
{
    scalar deltaT = fvMesh_.time().deltaT().value();
    const vectorField& oldPoints = fvMesh_.points();
    vectorField& motionPointUI = pointMotionU_.internalField();

    vectorField pointDisplacement(oldPoints.size(), vector::zero);

    forAll(motionPointUI, pointI)
    {
        motionPointUI[pointI] =
            motionPointUI[pointPoint_[pointI]]
           *(
                (
                    (motionDirection_ & oldPoints[pointI])
                  - (motionDirection_ & oldPoints[pointPoint2_[pointI]])
                )
               /(
                    (motionDirection_ & oldPoints[pointPoint_[pointI]])
                  - (motionDirection_ & oldPoints[pointPoint2_[pointI]])
                )
            );
    }

    tmp<pointField> tcurPoints
    (
        oldPoints + deltaT*pointMotionU_.internalField()
    );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}

// ************************************************************************* //
