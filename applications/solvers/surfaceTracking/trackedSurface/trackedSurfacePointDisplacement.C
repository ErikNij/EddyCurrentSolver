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

#include "trackedSurface.H"

#include "coordinateSystem.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


tmp<vectorField> trackedSurface::pointDisplacement(const scalarField& deltaH)
{
    const pointField& points = aMesh().patch().localPoints();
    const labelListList& pointFaces = aMesh().patch().pointFaces();

    controlPoints() += facesDisplacementDir()*deltaH;

// TEST: DEBUG | Additional debugging
    if (debug > 2)
    {
        mesh().write();

        writeVTKControlPoints();
    }

    if (correctCurvature_)
    {
        // Correct control points next to fixed patches
        forAll (fixedTrackedSurfacePatches_, patchI)
        {
            label fixedPatchID =
                aMesh().boundary().findPatchID
                (
                    fixedTrackedSurfacePatches_[patchI]
                );

            if (fixedPatchID == -1)
            {
                FatalErrorIn("trackedSurface::trackedSurface(...) : ")
                    << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                        << " defined in the trackedSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            const labelListList& fFaces = aMesh().patch().faceFaces();
            const vectorField& fCentres =
                aMesh().areaCentres().internalField();

            forAll (eFaces, edgeI)
            {
                const label& curFace = eFaces[edgeI];
                const labelList& curFaceFaces = fFaces[curFace];

                scalar H = 0.0;
                label counter = 0;

                forAll (curFaceFaces, faceI)
                {
                    label index = findIndex(eFaces, curFaceFaces[faceI]);

                    if (index == -1)
                    {
                        H +=
                            facesDisplacementDir()[curFaceFaces[faceI]]
                          & (
                                controlPoints()[curFaceFaces[faceI]]
                              - fCentres[curFaceFaces[faceI]]
                            );

                        counter++;
                    }
                }

                H /= counter;

                controlPoints()[curFace] =
                    fCentres[curFace]
                  + facesDisplacementDir()[curFace]*H;
            }
        }
    }


    tmp<vectorField> tdisplacement
    (
        new vectorField
        (
            points.size(),
            vector::zero
        )
    );

    vectorField& displacement = tdisplacement();


    // Calculate displacement of internal points
    const vectorField& pointNormals = aMesh().pointAreaNormals();
    const edgeList& edges = aMesh().patch().edges();
    labelList internalPoints = aMesh().internalPoints();

    forAll (internalPoints, pointI)
    {
        label curPoint = internalPoints[pointI];

        const labelList& curPointFaces = pointFaces[curPoint];

        vectorField lsPoints(curPointFaces.size(), vector::zero);

        for (label i=0; i<curPointFaces.size(); i++)
        {
            label curFace = curPointFaces[i];

            lsPoints[i] = controlPoints()[curFace];
        }

        vectorField pointAndNormal =
            lsPlanePointAndNormal
            (
                lsPoints,
                points[curPoint],
                pointNormals[curPoint]
            );

        vector& P = pointAndNormal[0];
        vector& N = pointAndNormal[1];

        displacement[curPoint] =
            pointsDisplacementDir()[curPoint]
           *((P - points[curPoint])&N)
           /(pointsDisplacementDir()[curPoint]&N);
    }


    // Mirror control points
    FieldField<Field, vector> patchMirrorPoints(aMesh().boundary().size());

// TEST: Why oldPoints?
//     // Old faMesh points
//     vectorField oldPoints(aMesh().nPoints(), vector::zero);
//     const labelList& meshPoints = aMesh().patch().meshPoints();
//     forAll (oldPoints, pI)
//     {
//         oldPoints[pI] =
//             mesh().oldPoints()[meshPoints[pI]];
//     }

    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         != processorFaPatch::typeName
        )
        {
            patchMirrorPoints.set
            (
                patchI,
                new vectorField
                (
                    aMesh().boundary()[patchI].faPatch::size(),
                    vector::zero
                )
            );

            vectorField N =
                aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

            vectorField Nr = N;

            const labelList& eFaces =
                aMesh().boundary()[patchI].edgeFaces();

            // Correct Nr according to specified contact angle
            if (contactAnglePtr_)
            {
                label ngbPolyPatchID =
                    aMesh().boundary()[patchI].ngbPolyPatchIndex();

                if (ngbPolyPatchID != -1)
                {
                    if
                    (
                        isA<wallFvPatch>(mesh().boundary()[ngbPolyPatchID])
                    )
                    {
                        scalarField& contactAngle =
                            contactAnglePtr_->boundaryField()[patchI];

                        scalarField rotAngle = 90 - contactAngle;

//                         rotAngle = average(rotAngle);

                        rotAngle *= M_PI/180.0;

                        vectorField rotationAxis(N.size(), vector::zero);

                        const vectorField& pEdgN =
                            aMesh().edgeAreaNormals().boundaryField()[patchI];

                        rotationAxis = (N^pEdgN);

                        const edgeList::subList patchEdges =
                            aMesh().boundary()[patchI].patchSlice(aMesh().edges());

                        forAll (rotationAxis, edgeI)
                        {
// TEST: Why oldPoints?
//                             vector e = patchEdges[edgeI].vec(oldPoints);
                            vector e = patchEdges[edgeI].vec(aMesh().points());

                            // Adjust direction
                            rotationAxis[edgeI] =
                                e*(e&rotationAxis[edgeI])
                            /mag((e&rotationAxis[edgeI]));
                        }
                        rotationAxis /= mag(rotationAxis) + SMALL;

                        vectorField rotationAxis2 = rotationAxis;
                        forAll (rotationAxis2, edgeI)
                        {
                            rotationAxis2[edgeI] =
                                (N[edgeI]^facesDisplacementDir()[eFaces[edgeI]]);

                            // Adjust direction
                            rotationAxis2[edgeI] =
                                rotationAxis2[edgeI]
                            *(rotationAxis2[edgeI]&rotationAxis[edgeI])
                            /mag((rotationAxis2[edgeI]&rotationAxis[edgeI]));
                        }
                        rotationAxis2 /= mag(rotationAxis2) + SMALL;

                        forAll (rotationAxis, edgeI)
                        {
                            vector NI = N[edgeI];
                            scalar rotAngleI = rotAngle[edgeI];
                            vector rotAxisI = rotationAxis[edgeI];

                            // Rodrigues' rotation formula
                            Nr[edgeI] = NI*cos(rotAngleI)
                            + rotAxisI*(rotAxisI & NI)*(1 - cos(rotAngleI))
                            + (rotAxisI^NI)*sin(rotAngleI);
                        }

                        Nr /= mag(Nr);

                        Nr = (rotationAxis^Nr);

                        Nr = (Nr^rotationAxis2);

                        Nr /= mag(Nr);
                    }
                }
            }

            const labelList peFaces =
                labelList::subList
                (
                    aMesh().edgeOwner(),
                    aMesh().boundary()[patchI].faPatch::size(),
                    aMesh().boundary()[patchI].start()
                );

            const labelList& pEdges = aMesh().boundary()[patchI];

            vectorField peCentres(pEdges.size(), vector::zero);
            forAll (peCentres, edgeI)
            {
                peCentres[edgeI] =
                    edges[pEdges[edgeI]].centre(points);
            }

            vectorField delta =
                vectorField(controlPoints(), peFaces)
            - peCentres;

//             patchMirrorPoints[patchI] =
//                 peCentres + ((I - 2*N*N)&delta);

            vectorField deltaN = N * (N & delta);

            vectorField deltaNr = Nr / (Nr&N) * mag(deltaN);

            patchMirrorPoints[patchI] =
                peCentres + delta + 2*deltaNr;

// TEST: DEBUG | Additional debugging
            if (debug > 2)
            {
                writeVTKpoints
                (
                    word("MirrorPoints")
                  + word(aMesh().boundary()[patchI].name()),
                    patchMirrorPoints[patchI]
                );
            }
        }
    }


    // Calculate displacement of boundary points
    labelList boundaryPoints = aMesh().boundaryPoints();

    const labelListList& edgeFaces = aMesh().patch().edgeFaces();
    const labelListList& pointEdges = aMesh().patch().pointEdges();

    forAll (boundaryPoints, pointI)
    {
        label curPoint = boundaryPoints[pointI];

        if (motionPointsMask()[curPoint] == 1)
        {
            // Calculating mirror points
            const labelList& curPointEdges = pointEdges[curPoint];

            vectorField mirrorPoints(2, vector::zero);

            label counter = -1;

            forAll (curPointEdges, edgeI)
            {
                label curEdge = curPointEdges[edgeI];

                if (edgeFaces[curEdge].size() == 1)
                {
                    label patchID = -1;
                    label edgeID = -1;
                    forAll (aMesh().boundary(), patchI)
                    {
                        const labelList& pEdges =
                            aMesh().boundary()[patchI];
                        label index = findIndex(pEdges, curEdge);
                        if (index != -1)
                        {
                            patchID = patchI;
                            edgeID = index;
                            break;
                        }
                    }

                    mirrorPoints[++counter] =
                        patchMirrorPoints[patchID][edgeID];
                }
            }

            // Calculating LS plane fit
            const labelList& curPointFaces = pointFaces[curPoint];

            vectorField lsPoints
            (
                curPointFaces.size() + mirrorPoints.size(),
                vector::zero
            );

            counter = -1;

            for (label i=0; i<curPointFaces.size(); i++)
            {
                label curFace = curPointFaces[i];

                lsPoints[++counter] = controlPoints()[curFace];
            }

            for (label i=0; i<mirrorPoints.size(); i++)
            {
                lsPoints[++counter] = mirrorPoints[i];
            }

            vectorField pointAndNormal =
                lsPlanePointAndNormal
                (
                    lsPoints,
                    points[curPoint],
                    pointNormals[curPoint]
                );

            vector& P = pointAndNormal[0];
            vector& N = pointAndNormal[1];

            displacement[curPoint] =
                pointsDisplacementDir()[curPoint]
               *((P - points[curPoint])&N)
               /(pointsDisplacementDir()[curPoint]&N);
        }
    }


    // Calculate displacement of axis point
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            forAll (wedgePatch.axisPoints(), pI)
            {
                label axisPoint = wedgePatch.axisPoints()[pI];

                displacement[axisPoint] =
                    pointsDisplacementDir()[axisPoint]
                   *(
                        pointsDisplacementDir()[axisPoint]
                      & (
                            controlPoints()[pointFaces[axisPoint][0]]
                          - points[axisPoint]
                        )
                    );
            }

            // ZT, 22-Sep-2014
//             if (wedgePatch.axisPoint() > -1)
//             {
//                 label axisPoint = wedgePatch.axisPoint();

//                 displacement[axisPoint] =
//                     pointsDisplacementDir()[axisPoint]
//                    *(
//                         pointsDisplacementDir()[axisPoint]
//                        &(
//                             controlPoints()[pointFaces[axisPoint][0]]
//                           - points[axisPoint]
//                         )
//                     );
//             }
        }
    }


    // Calculate displacement of processor patch points
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPointLabels =
                procPatch.pointLabels();

            FieldField<Field, vector> lsPoints(patchPointLabels.size());
            forAll (lsPoints, pointI)
            {
                lsPoints.set(pointI, new vectorField(0, vector::zero));
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll (nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                const labelList& curPointFaces = pointFaces[curPoint];

                lsPoints[curPatchPoint].setSize(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    lsPoints[curPatchPoint][faceI] = controlPoints()[curFace];
                }

#               include "boundaryProcessorFaPatchPoints.H"
            }

//             scalar lsPointsSize = 0;
//             forAll (lsPoints, pointI)
//             {
//                 lsPointsSize +=
//                     2*lsPoints[pointI].size()*sizeof(vector);
//             }

            // Parallel data exchange
            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
//                     lsPointsSize
                );

                toNeighbProc << lsPoints;
            }

            FieldField<Field, vector> ngbLsPoints(patchPointLabels.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
//                     lsPointsSize
                );

                fromNeighbProc >> ngbLsPoints;
            }

            forAll (nonGlobalPatchPoints, pointI)
            {
                label curPatchPoint =
                    nonGlobalPatchPoints[pointI];

                label curPoint =
                    patchPointLabels[curPatchPoint];

                label curNgbPoint = procPatch.neighbPoints()[curPatchPoint];

                vectorField allLsPoints
                (
                    lsPoints[curPatchPoint].size()
                  + ngbLsPoints[curNgbPoint].size(),
                    vector::zero
                );

                label counter = -1;
                forAll (lsPoints[curPatchPoint], pointI)
                {
                    allLsPoints[++counter] = lsPoints[curPatchPoint][pointI];
                }
                forAll (ngbLsPoints[curNgbPoint], pointI)
                {
                    allLsPoints[++counter] = ngbLsPoints[curNgbPoint][pointI];
                }

                vectorField pointAndNormal =
                    lsPlanePointAndNormal
                    (
                        allLsPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    );

                vector& P = pointAndNormal[0];
                vector& N = pointAndNormal[1];

                if (motionPointsMask()[curPoint] != 0)
                {
                    displacement[curPoint] =
                        pointsDisplacementDir()[curPoint]
                       *((P - points[curPoint])&N)
                       /(pointsDisplacementDir()[curPoint]&N);
                }
            }
        }
    }


    // Calculate displacement of global processor patch points
    if (aMesh().globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels =
            aMesh().globalData().sharedPointLabels();

        const labelList& addr = aMesh().globalData().sharedPointAddr();

        for (label k=0; k<aMesh().globalData().nGlobalPoints(); k++)
        {
            List<List<vector> > procLsPoints(Pstream::nProcs());

            label curSharedPointIndex = findIndex(addr, k);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelList& curPointFaces = pointFaces[curPoint];

                procLsPoints[Pstream::myProcNo()] =
                    List<vector>(curPointFaces.size());

                forAll (curPointFaces, faceI)
                {
                    label curFace = curPointFaces[faceI];

                    procLsPoints[Pstream::myProcNo()][faceI] =
                        controlPoints()[curFace];
                }
            }

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                label nAllPoints = 0;
                forAll (procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, vector::zero);

                label counter = 0;
                forAll (procLsPoints, procI)
                {
                    forAll (procLsPoints[procI], pointI)
                    {
                        allPoints[counter++] =
                            procLsPoints[procI][pointI];
                    }
                }

                vectorField pointAndNormal =
                    lsPlanePointAndNormal
                    (
                        allPoints,
                        points[curPoint],
                        pointNormals[curPoint]
                    );

                const vector& P = pointAndNormal[0];
                const vector& N = pointAndNormal[1];

                displacement[curPoint] =
                    pointsDisplacementDir()[curPoint]
                   *((P - points[curPoint])&N)
                   /(pointsDisplacementDir()[curPoint]&N);
            }
        }
    }

// TEST: DEBUG | Additional debugging
    if (debug > 2)
    {
        writeVTKpoints
        (
            word("displacement"),
            displacement
        );

        const pointField& oldPoints = mesh().allPoints();

        const labelList& meshPoints =
            mesh().boundaryMesh()[aPatchID()].meshPoints();

        vectorField newPoints(displacement.size(), vector::zero);

        forAll (newPoints, pointI)
        {
            newPoints[pointI] =
                oldPoints[meshPoints[pointI]]
              + displacement[pointI];
        }

        writeVTKpoints
        (
            word("newMeshPoints"),
            newPoints
        );
    }

    return tdisplacement;
}



tmp<vectorField> trackedSurface::lsPlanePointAndNormal
(
    const vectorField& points,
    const vector& origin,
    const vector& axis
) const
{
    // LS in local CS
    vector dir = (points[0] - origin);
    dir -= axis*(axis&dir);
    dir /= mag(dir);
    coordinateSystem cs("cs", origin, axis, dir);

    vectorField localPoints = cs.localPosition(points);
    scalarField W = 1.0/(mag(points - origin) + SMALL);

    scalarRectangularMatrix M
    (
        points.size(),
        3,
        0.0
    );

    for (label i=0; i<localPoints.size(); i++)
    {
        M[i][0] = localPoints[i].x();
        M[i][1] = localPoints[i].y();
        M[i][2] = 1.0;
    }

    scalarSquareMatrix MtM(3, 0.0);
    for (label i = 0; i < MtM.n(); i++)
    {
        for (label j = 0; j < MtM.m(); j++)
        {
            for (label k = 0; k < M.n(); k++)
            {
                MtM[i][j] += M[k][i]*M[k][j]*W[k];
            }
        }
    }

    scalarField MtR(3, 0);
    for (label i = 0; i < MtR.size(); i++)
    {
        for (label j = 0; j < M.n(); j++)
        {
            MtR[i] += M[j][i]*localPoints[j].z()*W[j];
        }
    }

    scalarSquareMatrix::LUsolve(MtM, MtR);

    vector n0 = vector(-MtR[0], -MtR[1], 1);
    n0 = cs.globalVector(n0);
    n0 /= mag(n0);

    vector p0 = vector(0, 0, MtR[2]);
    p0 = cs.globalPosition(p0);

    tmp<vectorField> pointAndNormal
    (
        new vectorField(2, vector::zero)
    );

    pointAndNormal()[0] = p0;
    pointAndNormal()[1] = n0;

    return pointAndNormal;
}


void trackedSurface::correctPointDisplacement
(
    const scalarField& sweptVolCorr,
    vectorField& displacement
)
{
    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const faceList& faces =
        aMesh().patch().localFaces();

    const pointField& points =
        aMesh().patch().localPoints();

    forAll (fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if (fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::correctPointDisplacement(...) : ")
                 << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        labelHashSet pointSet;

        forAll (eFaces, edgeI)
        {
            label curFace = eFaces[edgeI];

            const labelList& curPoints = faces[curFace];

            forAll (curPoints, pointI)
            {
                label curPoint = curPoints[pointI];
                label index = findIndex(pLabels, curPoint);

                if (index == -1)
                {
                    if (!pointSet.found(curPoint))
                    {
                        pointSet.insert(curPoint);
                    }
                }
            }
        }

        labelList corrPoints = pointSet.toc();

        labelListList corrPointFaces(corrPoints.size());

        forAll (corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;

            forAll (pFaces[curPoint], faceI)
            {
                label curFace = pFaces[curPoint][faceI];

                label index = findIndex(eFaces, curFace);

                if (index != -1)
                {
                    if (!faceSet.found(curFace))
                    {
                        faceSet.insert(curFace);
                    }
                }
            }

            corrPointFaces[pointI] = faceSet.toc();
        }

        forAll (corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll (curPointFaces, faceI)
            {
                const face& curFace = faces[curPointFaces[faceI]];

                label ptInFace = curFace.which(curPoint);
                label next = curFace.nextLabel(ptInFace);
                label prev = curFace.prevLabel(ptInFace);

                vector a = points[next] - points[curPoint];
                vector b = points[prev] - points[curPoint];
                const vector& c = pointsDisplacementDir()[curPoint];

                curDisp += 2*sweptVolCorr[curPointFaces[faceI]]/((a^b)&c);
            }

            curDisp /= curPointFaces.size();

            displacement[curPoint] =
                curDisp*pointsDisplacementDir()[curPoint];
        }
    }
}


// tmp<vectorField> trackedSurface::pointDisplacementID(const scalarField& deltaH)
// {
//     const pointField& points = aMesh().patch().localPoints();
//     const labelListList& pointFaces = aMesh().patch().pointFaces();
//
//     const vectorField& areaCentres = aMesh().areaCentres().internalField();
//
//     controlPoints() = areaCentres;
//     controlPoints() += facesDisplacementDir()*deltaH;
//
//     tmp<vectorField> tdisplacement
//     (
//         new vectorField
//         (
//             points.size(),
//             vector::zero
//         )
//     );
//
//     vectorField& displacement = tdisplacement();
//
//
//     forAll (pointFaces, pointI)
//     {
//         scalar weightsSum = 0.0;
//
//         const labelList& curPointFaces = pointFaces[pointI];
//
//         forAll (curPointFaces, faceI)
//         {
//             label curFace = curPointFaces[faceI];
//
//             scalar weight = 1.0/mag
//             (
//                 points[pointI]
//               - areaCentres[curFace]
//             );
//
//             displacement[pointI] += weight*areaCentres[curFace];
//
//             weightsSum += weight;
//         }
//
//         displacement[pointI] /= weightsSum;
//
//         displacement[pointI] -= points[pointI];
//     }
//
//
//     displacement = pointsDisplacementDir()
//         * (pointsDisplacementDir()&displacement);
//
//     forAll (motionPointsMask(), pointI)
//     {
//         displacement[pointI] *= motionPointsMask()[pointI];
//     }
//
//     return tdisplacement;
// }


// tmp<vectorField> trackedSurface::pointDisplacementSM()
// {
//     const pointField& points = aMesh().patch().localPoints();
//     const labelListList& pointFaces = aMesh().patch().pointFaces();


//     tmp<vectorField> tdisplacement
//     (
//         new vectorField
//         (
//             points.size(),
//             vector::zero
//         )
//     );

//     vectorField& displacement = tdisplacement();


//     forAll (pointFaces, pointI)
//     {
//         scalar weightsSum = 0.0;

//         const labelList& curPointFaces = pointFaces[pointI];

//         forAll (curPointFaces, faceI)
//         {
//             label curFace = curPointFaces[faceI];

//             scalar weight = 1.0/mag
//             (
//                 points[pointI]
//               - controlPoints()[curFace]
//             );

//             displacement[pointI] += weight*controlPoints()[curFace];

//             weightsSum += weight;
//         }

//         displacement[pointI] /= weightsSum;

//         displacement[pointI] -= points[pointI];
//     }


//     displacement = motionPointsMask()*
//         (pointsDisplacementDir()&displacement)*
//         pointsDisplacementDir();


//     return tdisplacement;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
