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

#include "liquidFilmParallelFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidFilmParallelFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        liquidFilmParallelFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidFilmParallelFvMotionSolver::liquidFilmParallelFvMotionSolver
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
    motionDirection_(lookup("motionDirection")),
    freeSurfacePatchID_(mesh.boundaryMesh().findPatchID("freeSurface")),
    freeSurfaceZoneID_(mesh.faceZones().findZoneID("freeSurfaceZone")),
    procToGlobalFZmap_()
{
    motionDirection_ /= mag(motionDirection_) + SMALL;

    boundBox box(mesh.points());

    max_ = (motionDirection_&box.max());
    min_ = (motionDirection_&box.min());

    // Check free surface patch
    if(freeSurfacePatchID_ == -1)
    {
        FatalErrorIn
        (
            "liquidFilmParallelFvMotionSolver(...)"
        )
            << "Free surface patch not defined." << abort(FatalError);
    }

    // Check free surface face zone
    if(freeSurfaceZoneID_ == -1)
    {
        FatalErrorIn
        (
            "liquidFilmParallelFvMotionSolver(...)"
        )
            << "Free surface face zone not defined." << abort(FatalError);
    }

    const vectorField& points = fvMesh_.points();

    vectorField fzLocalPoints =
        mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

    forAll(pointPoint_, pointI)
    {
        scalar minDist = GREAT;

        forAll(fzLocalPoints, lpI)
        {
            vector R = fzLocalPoints[lpI] - points[pointI];
            R -= motionDirection_*(motionDirection_&R);
            scalar curDist = mag(R);

            if (curDist < minDist)
            {
                minDist = curDist;
                pointPoint_[pointI] = lpI;
            }
        }
    }

    // face zona map

    procToGlobalFZmap_.setSize(fzLocalPoints.size(), -1);

    if(Pstream::parRun())
    {
        vectorField fzGlobalPoints =
            mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

        //- set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints *= 0.0;
        }

        //- pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        //- now every proc has the master's list of FZ points
       	//- every proc must now find the mapping from their local FZ points to
        //- the global FZ points

        vectorField fzLocalPoints =
            mesh.faceZones()[freeSurfaceZoneID_]().localPoints();

        forAll(fzGlobalPoints, globalPointI)
        {
            forAll(fzLocalPoints, procPointI)
            {
                if
                (
                    mag
                    (
                        fzLocalPoints[procPointI]
                      - fzGlobalPoints[globalPointI]
                    )
                  < SMALL
                )
                {
                    procToGlobalFZmap_[globalPointI] = procPointI;
                    break;
                }
            }
        }

        forAll(procToGlobalFZmap_, globalPointI)
        {
            if (procToGlobalFZmap_[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "liquidFilmParallelFvMotionSolver(...)"
                )
                    << "Global FZ map is not correct"
                        << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liquidFilmParallelFvMotionSolver::~liquidFilmParallelFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::liquidFilmParallelFvMotionSolver::curPoints() const
{
    scalar deltaT = fvMesh_.time().deltaT().value();
    const vectorField& oldPoints = fvMesh_.points();
    vectorField& motionPointUI = pointMotionU_.internalField();

    const labelList fzMeshPoints =
        fvMesh_.faceZones()[freeSurfaceZoneID_]().meshPoints();

    // Displacement for face zone points
    vectorField fzMotionPointU(fzMeshPoints.size(), vector::zero);

    //- inter-proc points are shared by multiple procs
    //- pointNumProc is the number of procs which a point lies on
    scalarField pointNumProcs(fzMeshPoints.size(), 0);

    forAll(fzMotionPointU, globalPointI)
    {
        label localPoint = procToGlobalFZmap_[globalPointI];

        if
        (
            fzMeshPoints[localPoint] < fvMesh_.nPoints()
//             localPoint
//           < fvMesh_.boundaryMesh()[freeSurfacePatchID_].localPoints().size()
        )
        {
            label procPoint = fzMeshPoints[localPoint];

            fzMotionPointU[globalPointI] = motionPointUI[procPoint];

            pointNumProcs[globalPointI] = 1;
        }
    }

    reduce(fzMotionPointU, sumOp<vectorField>());
    reduce(pointNumProcs, sumOp<scalarField>());

    Pout << "pointNumProcs" << min(pointNumProcs) << ", "
        << max(pointNumProcs) << endl;

    //- now average the motionPointU between all procs
    fzMotionPointU /= pointNumProcs;

    Pout << fzMotionPointU[0] << ", "
        << fzMotionPointU[fzMotionPointU.size()-1] << endl;

    //- the globalFZnewPoints now contains the correct FZ displacement in
    //- a global order, now convert them back into the local proc order

    vectorField fzLocalMotionPointU(fzMeshPoints.size(), vector::zero);

    forAll(fzMotionPointU, globalPointI)
    {
        label localPoint = procToGlobalFZmap_[globalPointI];

        fzLocalMotionPointU[localPoint] = fzMotionPointU[globalPointI];
    }

    const vectorField fzLocalPoints =
        fvMesh_.faceZones()[freeSurfaceZoneID_]().localPoints();

    forAll(motionPointUI, pointI)
    {
        if
        (
            (motionDirection_ & oldPoints[pointI])
         <= (motionDirection_ & fzLocalPoints[pointPoint_[pointI]])
        )
        {
            motionPointUI[pointI] =
                fzLocalMotionPointU[pointPoint_[pointI]]
               *(
                    (
                        (motionDirection_ & oldPoints[pointI])
                      - min_
                    )
                   /(
                        (
                            motionDirection_
                          & fzLocalPoints[pointPoint_[pointI]]
                        )
                      - min_
                    )
                );
        }
        else
        {
            motionPointUI[pointI] =
                fzLocalMotionPointU[pointPoint_[pointI]]
               *(
                    (
                        max_
                      - (motionDirection_ & oldPoints[pointI])
                    )
                   /(
                        max_
                      - (
                            motionDirection_
                          & fzLocalPoints[pointPoint_[pointI]]
                        )
                    )
                );
        }
    }

    vectorField newPoints = fvMesh_.points();

    forAll(motionPointUI, pointI)
    {
        newPoints[pointI] += deltaT*motionPointUI[pointI];
    }

    Pout << "End" << endl;
    ::sleep(10);

//     forAll(fzLocalMotionPointU, pointI)
//     {
//         if (fzMeshPoints[pointI] >= fvMesh_.nPoints())
//         {
//             newPoints[fzMeshPoints[pointI]] +=
//                 deltaT*fzLocalMotionPointU[pointI];
//         }
//     }


//     Pout << "End2" << endl;
//     ::sleep(10);

    tmp<pointField> tcurPoints(newPoints);

//     tmp<pointField> tcurPoints
//     (
//         oldPoints + deltaT*pointMotionU_.internalField()
//     );

    twoDCorrectPoints(tcurPoints());

//     Pout << "End3" << endl;
//     ::sleep(10);

//     Pout << tcurPoints() << endl;

    return tcurPoints;
}

// ************************************************************************* //
