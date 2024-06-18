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

#include "liquidFilmFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidFilmFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        liquidFilmFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidFilmFvMotionSolver::liquidFilmFvMotionSolver
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
    max_(0),
    min_(0),
    motionDirection_(lookup("motionDirection"))
{
    motionDirection_ /= mag(motionDirection_) + SMALL;

    boundBox box(mesh.points());

    max_ = (motionDirection_&box.max());
    min_ = (motionDirection_&box.min());

    // Detect the free surface patch

    label aPatchID = -1;

    forAll (fvMesh_.boundary(), patchI)
    {
        if(fvMesh_.boundary()[patchI].name() == "freeSurface")
        {
            aPatchID = patchI;

            Info<< "Found free surface patch. ID: " << aPatchID
                << endl;
        }
    }

    if(aPatchID == -1)
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }

    label bPatchID = -1;

    forAll (fvMesh_.boundary(), patchI)
    {
        if(fvMesh_.boundary()[patchI].name() == "freeSurfaceShadow")
        {
            bPatchID = patchI;

            Info<< "Found free surface shadow patch. ID: " << bPatchID
                << endl;
        }
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidFilmFvMotionSolver::~liquidFilmFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::liquidFilmFvMotionSolver::curPoints() const
{
    scalar deltaT = fvMesh_.time().deltaT().value();
    const vectorField& oldPoints = fvMesh_.points();
    vectorField& motionPointUI = pointMotionU_.internalField();

    vectorField pointDisplacement(oldPoints.size(), vector::zero);

    forAll(motionPointUI, pointI)
    {
        if
        (
            (motionDirection_ & oldPoints[pointI])
         <= (motionDirection_ & oldPoints[pointPoint_[pointI]])
        )
        {
            motionPointUI[pointI] =
                motionPointUI[pointPoint_[pointI]]
               *(
                    (
                        (motionDirection_ & oldPoints[pointI])
                      - min_
                    )
                   /(
                        (motionDirection_ & oldPoints[pointPoint_[pointI]])
                      - min_
                    )
                );
        }
        else
        {
            motionPointUI[pointI] =
                motionPointUI[pointPoint_[pointI]]
               *(
                    (
                        max_
                      - (motionDirection_ & oldPoints[pointI])
                    )
                   /(
                        max_
                      - (motionDirection_ & oldPoints[pointPoint_[pointI]])
                    )
                );
        }
    }

    tmp<pointField> tcurPoints
    (
        oldPoints + deltaT*pointMotionU_.internalField()
    );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}

// ************************************************************************* //
