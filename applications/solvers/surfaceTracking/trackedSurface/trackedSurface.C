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

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "twoDPointCorrector.H"

#include "tetFemMatrices.H"
#include "faceTetPolyPatch.H"
#include "tetPointFields.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"

#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "fixedGradientFaPatchFields.H"

// TODO
#include "coordinateSystem.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(trackedSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void trackedSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
// TEST: Move always from start
    deleteDemandDrivenData(points0Ptr_);
// TEST: Move always from start
    deleteDemandDrivenData(total0DisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
// TEST: Sub-mesh
    deleteDemandDrivenData(aSubMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
    deleteDemandDrivenData(muEffFluidAvalPtr_);
    deleteDemandDrivenData(muEffFluidBvalPtr_);
    deleteDemandDrivenData(contactAnglePtr_);
    deleteDemandDrivenData(temperaturePtr_);
    deleteDemandDrivenData(concentrationPtr_);
    deleteDemandDrivenData(surfaceTensionForcePtr_);
    deleteDemandDrivenData(nGradUnPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void trackedSurface::initProperties()
{
    updateProperties();
}


void trackedSurface::initCheckPointNormalsCorrection()
{
    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll (pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if (patchID == -1)
        {
            FatalErrorIn
            (
                "trackedSurface::initCheckPointNormalsCorrection(...) : "
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Check correction for patches with specified contact angle
    if (contactAnglePtr_)
    {
        forAll (aMesh().boundary(), patchI)
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
                    if (correction[patchI])
                    {
                        FatalErrorIn
                        (
                            "trackedSurface::initCheckPointNormalsCorrection(...) : "
                        )   << "Point normal correction is activated for patch "
                            << aMesh().boundary()[patchI].name() << ". "
                            << "This is considered fatal as a contact angle is "
                            << "beeing specified here, too."
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


void trackedSurface::initCheckSurfacePatches()
{
    // Detect the surface patch
    forAll (mesh().boundary(), patchI)
    {
        if (mesh().boundary()[patchI].name() == prefix_)
        {
            aPatchID_ = patchI;

            if (debug)
            {
                Info << "trackedSurface::initCheckSurfacePatches() : "
                    << "Found surface patch. ID: " << aPatchID_
                    << endl;
            }
        }
    }

    if (aPatchID() == -1)
    {
        FatalErrorIn("trackedSurface::initCheckSurfacePatches(...) : ")
            << "Surface patch not defined.  Please make sure that "
                << " the surface patches is named as trackedSurface"
                << abort(FatalError);
    }


    // Detect the surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if (mesh().boundary()[patchI].name() == prefix_ + "Shadow")
            {
                bPatchID_ = patchI;

                if (debug)
                {
                    Info << "trackedSurface::initCheckSurfacePatches() : "
                        << "Found surface shadow patch. ID: " << bPatchID_
                        << endl;
                }
            }
        }

        if (bPatchID() == -1)
        {
            FatalErrorIn("trackedSurface::initCheckSurfacePatches(...) : ")
                << "Surface shadow patch not defined. "
                    << "Please make sure that the surface shadow patch "
                    << "is named as trackedSurfaceShadow."
                    << abort(FatalError);
        }
    }
}


void trackedSurface::initMotionPointMask()
{
    // Mark surface boundary points
    // which belonge to processor patches
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll (patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }

    // Mark fixed surface boundary points
    forAll (fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if (fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::initMotionPointMask(...) : ")
                << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll (patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }

    // Mark surface boundary point
    // at the axis of 2-D axisymmetic cases
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
                motionPointsMask()[wedgePatch.axisPoints()[pI]] = 0;
            }

//             if (wedgePatch.axisPoint() > -1)
//             {
//                 motionPointsMask()[wedgePatch.axisPoint()] = 0;

//                 Info << "Axis point: "
//                     << wedgePatch.axisPoint()
//                     << " vector: "
//                     << aMesh().points()[wedgePatch.axisPoint()] << endl;
//             }
        }
    }
}


void trackedSurface::initControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    vectorField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();

    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll (faces, faceI)
    {
        face curFace = faces[faceI];
        sweptVol[faceI] = -curFace.sweptVol(points, newPoints);
    }

    vectorField faceAreaNormals(faces.size(), vector::zero);

    forAll (faceAreaNormals, faceI)
    {
        face curFace = faces[faceI];
        faceAreaNormals[faceI] = curFace.normal(newPoints);
    }

    forAll (deltaH, faceI)
    {
        scalar AnI = faceAreaNormals[faceI] & facesDisplacementDir()[faceI];

        deltaH[faceI] = sweptVol[faceI]/ (AnI + SMALL);

        if (AnI < SMALL)
        {
            FatalErrorIn("trackedSurface::initControlPointsPosition(...)")
                << "Something is probably wrong with the specified motion direction"
                    << abort(FatalError);

        }
    }

    calcAddDeltaHcorrection(deltaH);

    displacement = pointDisplacement(deltaH);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void trackedSurface::resetControlPoints(bool force)
{
    if (force || resetControlPoints_)
    {
        controlPoints() = aMesh().areaCentres().internalField();
    }
}


void trackedSurface::adjustDisplacementDirections()
{
    if (normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correct point displacement direction
        // at the "axis" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll (aMesh().boundary(), patchI)
        {
            if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

                vector axis = wedgePatch.axis();

                label centerLinePatchID =
                    aMesh().boundary().findPatchID("axis");

                if (centerLinePatchID != -1)
                {
                    const labelList& pointLabels =
                        aMesh().boundary()[centerLinePatchID].pointLabels();

                    forAll (pointLabels, pointI)
                    {
                        vector dir =
                            pointsDisplacementDir()[pointLabels[pointI]];

                        dir = (dir&axis)*axis;
                        dir /= mag(dir);

                        pointsDisplacementDir()[pointLabels[pointI]] = dir;
                    }
                }
                else
                {
                    WarningIn("trackedSurface::adjustDisplacementDirections()")
                        << "Centerline polyPatch does not exist. "
                        << "Surface points displacement directions "
                        << "will not be corrected at the axis (axis)"
                        << endl;
                }

                break;
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points postion
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TEST: Volume conservation
// TODO: Check this implementation in parallel
void trackedSurface::correctInterfacePhi(scalarField& interfacePhi)
{
    if(correctVolume_)
    {
        if (debug)
        {
            Info << "trackedSurface::correctInterfacePhi() : "
                << "Correct interface phi to account"
                    << " for domain volume changes."
                    << endl;
        }

        scalar phiAsumPatchesNotAB = 0.0;

        forAll (phi().boundaryField(), patchI)
        {
            const fvsPatchScalarField& phip = phi().boundaryField()[patchI];
            const fvPatchScalarField& fip = fluidIndicator().boundaryField()[patchI];

            if (patchI != aPatchID())
            {
                if
                (
                    (twoFluids() && patchI != bPatchID())
                && !twoFluids()
                )
                {
                    forAll (phip, faceI)
                    {
                        phiAsumPatchesNotAB += fip[faceI] * phip[faceI];
                    }
                }
            }
        }

        scalar vol =  fvc::domainIntegrate(fluidIndicator()).value();

        scalar phiAsumDomain =
            (vol - vol0_) / DB().deltaT().value();

// TODO: Think about best weights for distributing the volume correction.
        const scalarField& weights = aMesh().S();
//         const scalarField weights = mag(interfacePhi);

        scalarField domainPhi
        (
            weights * (phiAsumDomain - phiAsumPatchesNotAB) / gSum(weights)
        );

        scalar volDiffDomainRel =
            (phiAsumDomain - phiAsumPatchesNotAB)
        * DB().deltaT().value()
        / vol0_;

        if (debug > 1)
        {
            Info << "trackedSurface::correctInterfacePhi() : "
                << "Statistics:" << endl
                    << "  phiAsumPatchesNotAB = " << phiAsumPatchesNotAB << endl
                    << "  phiAsumDomain = " << phiAsumDomain << endl
                    << "  volDiffDomainRel = " << 100*volDiffDomainRel << " %"
                    << endl;
        }


// TODO: How to inject domainPhi best? Just adding it to meshPhi
//       dramatically decreases the convergence behaviour. Thus, in
//       the current implementation, phi is corrected, directly.
        interfacePhi -= domainPhi;
    }
}

tmp<scalarField> trackedSurface::calcSweptVolCorr(const scalarField& interfacePhi)
{
    tmp<scalarField> tsweptVolCorr
    (
        new scalarField
        (
            interfacePhi
          - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()]
        )
    );

    scalarField& sweptVolCorr = tsweptVolCorr();

    word ddtScheme
    (
        mesh().schemesDict().ddtScheme
        (
            "ddt(" + rho().name() + ',' + U().name()+')'
        )
    );

    if
    (
        ddtScheme
     == fv::CrankNicolsonDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= (1.0/2.0)*DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        sweptVolCorr *= DB().deltaT().value();
    }
    else if
    (
        ddtScheme
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        if (DB().timeIndex() == 1)
        {
            sweptVolCorr *= DB().deltaT().value();
        }
        else
        {
            sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
        }
    }
    else
    {
        FatalErrorIn("trackedSurface::calcSweptVolCorr()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    return tsweptVolCorr;
}


tmp<scalarField> trackedSurface::calcDeltaH(const scalarField& sweptVolCorr)
{
    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    tmp<scalarField> tdeltaH
    (
        new scalarField
        (
            sweptVolCorr.size(),
            0
        )
    );

    scalarField& deltaH = tdeltaH();

    deltaH = sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    calcAddDeltaHcorrection(deltaH);

    return tdeltaH;
}


void trackedSurface::calcAddDeltaHcorrection(scalarField& deltaH)
{
    forAll (fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if (fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::calcDeltaH(...) : ")
                << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll (eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }
}


tmp<vectorField> trackedSurface::calcDisplacement(const scalarField& interfacePhi)
{
    scalarField sweptVolCorr = calcSweptVolCorr(interfacePhi);
    scalarField deltaH = calcDeltaH(sweptVolCorr);

    tmp<vectorField> tdisplacement = pointDisplacement(deltaH);

    vectorField& displacement = tdisplacement();

    if (correctDisplacement_)
    {
        correctPointDisplacement(sweptVolCorr, displacement);
    }

    return tdisplacement;
}


tmp<pointField> trackedSurface::calcNewMeshPoints(const vectorField& displacement)
{
    tmp<pointField> tnewMeshPoints
    (
        new pointField(mesh().allPoints())
    );

    pointField& newMeshPoints = tnewMeshPoints();

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if (twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        vectorField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI];
        }

        Pout << gMax(mag(displacement)) << ", "
            << gMax(mag(displacementB)) << endl;
    }

    twoDPointCorrector twoDPointCorr(mesh());
    twoDPointCorr.correctPoints(newMeshPoints);

    return tnewMeshPoints;
}


tmp<pointField> trackedSurface::calcNewMeshPoints(const scalarField& interfacePhi)
{
    tmp<vectorField> tdisplacement = calcDisplacement(interfacePhi);

    vectorField& displacement = tdisplacement();

    tmp<pointField> tnewMeshPoints = calcNewMeshPoints(displacement);

    tdisplacement.clear();

    return tnewMeshPoints;
}


tmp<pointField> trackedSurface::calcNewMeshPoints()
{
    tmp<pointField> tnewMeshPoints
    (
        new pointField(mesh().allPoints())
    );

    pointField& newMeshPoints = tnewMeshPoints();

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (totalDisplacement(), pointI)
    {
        newMeshPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI];
    }

    twoDPointCorrector twoDPointCorr(mesh());
    twoDPointCorr.correctPoints(newMeshPoints);

    return tnewMeshPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void trackedSurface::movePoints(const scalarField& interfacePhi)
{
    vectorField displacement = calcDisplacement(interfacePhi);

    totalDisplacement() += displacement;

    if (total0Update_)
    {
        total0Displacement() += displacement;
    }

    pointField newMeshPoints = calcNewMeshPoints(displacement);

    mesh().movePoints(newMeshPoints);

    moveUpdate();
}


void trackedSurface::moveMeshPoints(const pointField& displacement)
{
    bool feMotionSolver =
        mesh().objectRegistry::foundObject<tetPointVectorField>
        (
            "motionU"
        );

    bool fvMotionSolver =
        mesh().objectRegistry::foundObject<pointVectorField>
        (
            "pointMotionU"
        );

    if (feMotionSolver)
    {
        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                displacement/DB().deltaT().value()
            );

        if (twoFluids_)
        {
            vectorField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );

            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch ==
                tppiBPatch.pointToPointInterpolate
                (
                    displacement/DB().deltaT().value()
                );
        }

        motionU.correctBoundaryConditions();
    }
    else if (fvMotionSolver)
    {
        pointVectorField& pointMotionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& pointMotionUaPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                pointMotionU.boundaryField()[aPatchID()]
            );

        pointMotionUaPatch ==
            displacement/DB().deltaT().value();

        if (twoFluids_)
        {
            vectorField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );

            fixedValuePointPatchVectorField& pointMotionUbPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    pointMotionU.boundaryField()[bPatchID()]
                );

            pointMotionUbPatch ==
                displacementB/DB().deltaT().value();
        }

        pointMotionU.correctBoundaryConditions();
    }

    mesh().update();
}


void trackedSurface::moveMeshPoints()
{
    pointField resetMeshPoints = calcNewMeshPoints();

// TEST: Move always from start
    if (total0Update_)
    {
        // Reset points to 0 position
        mesh().movePoints(points0());

        // Update points according to new surface shape
        moveMeshPoints(total0Displacement());

        // Store new point positions
        pointField newMeshPoints = mesh().allPoints();

        // Reset points to position of last time step
        mesh().movePoints(resetMeshPoints);

        // Move points back to new position
        mesh().movePoints(newMeshPoints);
    }
    else
    {
        // Reset points to position of last time step
        mesh().movePoints(resetMeshPoints);

        // Update points according to new surface shape
        moveMeshPoints(totalDisplacement());
    }

    deleteDemandDrivenData(totalDisplacementPtr_);
}


// TEST: Sub-mesh
void trackedSurface::moveFaSubMesh()
{
    if (aSubMeshPtr_)
    {
        if (debug)
        {
            Info << "trackedSurface::moveFaSubMesh() : "
                << "Moving/Updating triangulated sub-mesh."
                    << endl;
        }

        aSubMesh().movePoints();
    }
}


template<class Type>
void
trackedSurface::moveCorrectedPatchSubMesh
(
    GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    forAll (vf.boundaryField(), patchI)
    {
        if
        (
            (
                vf.boundaryField()[patchI].type()
                == fixedGradientCorrectedFvPatchField<Type>::typeName
            )
        ||
            (
                vf.boundaryField()[patchI].type()
                == fixedValueCorrectedFvPatchField<Type>::typeName
            )
        ||
            (
                vf.boundaryField()[patchI].type()
                == zeroGradientCorrectedFvPatchField<Type>::typeName
            )
        )
        {
            correctedFvPatchField<Type>& vfA =
                refCast<correctedFvPatchField<Type> >
                (
                    vf.boundaryField()[patchI]
                );

            vfA.movePatchSubMesh();
        }
    }
}


void trackedSurface::moveCorrectedPatchSubMeshes()
{
    moveCorrectedPatchSubMesh(U());

    moveCorrectedPatchSubMesh(p());

    if (TPtr_)
    {
        moveCorrectedPatchSubMesh(T());
    }
}


void trackedSurface::moveUpdate()
{
    // Move faSubMesh
    moveFaSubMesh();

    // Move correctedFvPatchFields
    moveCorrectedPatchSubMeshes();

// TODO
    correctContactLinePointNormals();

// TODO
    if (correctPointNormals_)
    {
        correctPointNormals();
    }

// TODO
    if (correctCurvature_)
    {
        correctCurvature();
    }

// TODO
    if (smoothCurvature_)
    {
        smoothCurvature();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

trackedSurface::trackedSurface
(
    dynamicFvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub,
    volScalarField& Pb,
    surfaceScalarField& sfPhi,
    volScalarField* TbPtr,
    volScalarField* cbPtr,
    const uniformDimensionedVectorField* gPtr,
    const twoPhaseMixture* transportPtr,
    const incompressible::turbulenceModel* turbulencePtr,
    const volScalarField* p0Ptr,
    word prefix
)
:
    IOdictionary
    (
        IOobject
        (
            prefix + "Properties",
            m.time().constant(),
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    prefix_(prefix),
    Prefix_(word(toupper(prefix_[0])) + word(prefix_.substr(1))),
// TEST: Volume conservation
    vol0_(0.0),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    phi_(sfPhi),
    TPtr_(TbPtr),
    cPtr_(cbPtr),
    gPtr_(gPtr),
    transportPtr_(transportPtr),
    turbulencePtr_(turbulencePtr),
    p0Ptr_(p0Ptr),
    curTimeIndex_(m.time().timeIndex()),
    meshMoved_(false),
    twoFluids_(false),
    normalMotionDir_(true),
    motionDir_(vector::zero),
    cleanInterface_(true),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_(word(), dimMass/dimLength/dimTime, 0),
    muFluidB_(word(), dimMass/dimLength/dimTime, 0),
    rhoFluidA_(word(), dimDensity, 0),
    rhoFluidB_(word(), dimDensity, 0),
    kFluidA_(word(), dimThermalConductivity, 0.0),
    kFluidB_(word(), dimThermalConductivity, 0.0),
    CpFluidA_(word(), dimSpecificHeatCapacity, 0.0),
    CpFluidB_(word(), dimSpecificHeatCapacity, 0.0),
    g_(word(), dimLength/pow(dimTime,2), vector::zero),
    cleanInterfaceSurfTension_(word(), dimMass/pow(dimTime,2), 0),
    fixedTrackedSurfacePatches_(wordList(0)),
    pointNormalsCorrectionPatches_(wordList(0)),
    nTrackedSurfCorr_(1),
// TEST: Fixed interface
    fixedInterface_(false),
// TEST: Sub-mesh
    useSubMesh_(false),
    correctVolume_(false),
// TEST: Move always from start
    total0Update_(false),
    resetControlPoints_(false),
    smoothing_(false),
    freeContactAngle_(false),
    correctPointNormals_(false),
    correctDisplacement_(false),
    correctCurvature_(false),
    smoothCurvature_(false),
    curvExtrapOrder_(0),
    fvcNGradUn_(false),
    implicitCoupling_(false),
    interfaceDeformationLimit_(0),
    interpolatorABPtr_(NULL),
    interpolatorBAPtr_(NULL),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    totalDisplacementPtr_(NULL),
// TEST: Move always from start
    points0Ptr_(NULL),
// TEST: Move always from start
    total0DisplacementPtr_(NULL),
    aMeshPtr_(NULL),
// TEST: Sub-mesh
    aSubMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    muEffFluidAvalPtr_(NULL),
    muEffFluidBvalPtr_(NULL),
    contactAnglePtr_(NULL),
    temperaturePtr_(NULL),
    concentrationPtr_(NULL),
    surfaceTensionForcePtr_(NULL),
    nGradUnPtr_(NULL)
{
    if (debug)
    {
        Info << "Surface prefix is: " << prefix_ << endl;
    }

    if (Switch(this->lookupOrDefault("fixedInterface", false)))
    {
        fixedInterface_ = true;

        twoFluids_ = false;
        normalMotionDir_ = true;
        fixedTrackedSurfacePatches_ = wordList(0);
        pointNormalsCorrectionPatches_ = wordList(0);
        freeContactAngle_ = true;

    }
    else
    {
        twoFluids_ = Switch(this->lookup("twoFluids"));

        normalMotionDir_ = Switch(this->lookup("normalMotionDir"));

        fixedTrackedSurfacePatches_ =
            wordList(this->lookup("fixed" + Prefix_ + "Patches"));

        pointNormalsCorrectionPatches_ =
            wordList(this->lookup("pointNormalsCorrectionPatches"));
    }

    // Init properties
    initProperties();

// TEST: Volume conservation
    vol0_ = fvc::domainIntegrate(fluidIndicator()).value();

    // Init motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(this->lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Make contact angle if necessary
    if
    (
        IOobject
        (
            "contactAngle",
            DB().timeName(),
            mesh(),
            IOobject::MUST_READ
        ).headerOk()
     || freeContactAngle_
    )
    {
// TODO
        makeContactAngle();
        updateContactAngle();
    }

    initCheckPointNormalsCorrection();

    initCheckSurfacePatches();

    initMotionPointMask();

    if (MarangoniStress())
    {
// TODO: Check with cPtr_
        // Check for multiple Marangoni effects
        if (!cleanInterface_ && TPtr_)
        {
            FatalErrorIn("trackedSurface::trackedSurface(...) : ")
                << "Marangoni effect due to both "
                    << "surfactant concentration gradient "
                    << "and temperature gradient is not implemented"
                    << abort(FatalError);
        }

        // Check for multiple Marangoni effects
        if (!cleanInterface_ && cPtr_)
        {
            FatalErrorIn("trackedSurface::trackedSurface(...) : ")
                << "Marangoni effect due to both "
                    << "surfactant concentration gradient "
                    << "and concentration gradient is not implemented"
                    << abort(FatalError);
        }

        // Make surfactant
        if (!cleanInterface_)
        {
            Info << "Surfactant-driven Marangoni effect enabled" << endl;
            surfactant();
            surfactantConcentration();
        }

        // Make surface temperature field
        if (TPtr_)
        {
            Info << "Temperature-driven Marangoni effect enabled" << endl;
            temperature();
        }

        // Make surface concentration field
        if (cPtr_)
        {
            Info << "Concentration-driven Marangoni effect enabled" << endl;
            concentration();
        }
    }

    if (!fixedInterface_)
    {
        // Init total displacement
        if
        (
            IOobject
            (
                "totalDisplacement",
                DB().timeName(),
                mesh(),
                IOobject::MUST_READ
            ).headerOk()
        )
        {
// TODO
            makeTotalDisplacement();
        }

// TEST: Move always from start
        if (total0Update_)
        {
            points0Ptr_ = new pointField(m.allPoints());

// TODO
//             readTotal0Displacement();
        }

// TODO
        // Init control points position
        initControlPointsPosition();

        // Clear geometry
        aMesh().movePoints();

// TEST: Sub-mesh
        if (useSubMesh_)
        {
            if (debug)
            {
                Info << "trackedSurface::trackedSurface(...) : "
                    << "Creating and activating triangulated sub-mesh."
                        << endl;
            }

            makeFaSubMesh();
            aSubMesh().movePoints();
        }

        // Contact angle correction
        correctContactLinePointNormals();
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

trackedSurface::~trackedSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// TODO: Why not use .lookupOrDefault() ????
void trackedSurface::updateProperties()
{
    if (transportPtr_)
    {

        rhoFluidA_ = transport().rho1();
        muFluidA_ = rhoFluidA_ * dimensionedScalar
        (
            transport().nuModel1().viscosityProperties().lookup("nu")
        );

        rhoFluidB_ = transport().rho2();
        muFluidB_ = rhoFluidB_ * dimensionedScalar
        (
            transport().nuModel2().viscosityProperties().lookup("nu")
        );

        cleanInterfaceSurfTension_ =
            dimensionedScalar(transport().lookup("sigma"));
    }
    else
    {
        rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));
        muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));

        rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));
        muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));

        cleanInterfaceSurfTension_ =
            dimensionedScalar(this->lookup("surfaceTension"));
    }

    if (TPtr_)
    {
        kFluidA_ = dimensionedScalar(this->lookup("kFluidA"));
        kFluidB_ = dimensionedScalar(this->lookup("kFluidB"));
        CpFluidA_ = dimensionedScalar(this->lookup("CpFluidA"));
        CpFluidB_ = dimensionedScalar(this->lookup("CpFluidB"));
    }

    if (gPtr_)
    {
        g_ = dimensionedVector(*gPtr_);
    }
    else
    {
        g_ = dimensionedVector(this->lookup("g"));
    }

    // Read surfactant setting
    cleanInterface_ =
        Switch(this->lookupOrDefault("cleanInterface", true));

    // Read correctors count
    nTrackedSurfCorr_ =
        int(this->lookupOrDefault("n" + Prefix_ + "Correctors", 1));

    // Check if sub-mesh switch is set
    if (this->found("useSubMesh"))
    {
        useSubMesh_ = Switch(this->lookup("useSubMesh"));
    };

    // Check if volume correction switch is set
    if (this->found("correctVolume"))
    {
        correctVolume_ = Switch(this->lookup("correctVolume"));
    };

    // Check if total0Update switch is set
    if (this->found("total0Update"))
    {
        total0Update_ = Switch(this->lookup("total0Update"));
    };

    // Check if reset control points switch is set
    if (this->found("resetControlPoints"))
    {
        resetControlPoints_ = Switch(this->lookup("resetControlPoints"));
    };

    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
    };

    // Check if freeContactAngle switch is set
    if (fixedInterface_)
    {
        freeContactAngle_ = true;
    }
    else
    {
        if (this->found("freeContactAngle"))
        {
            freeContactAngle_ = Switch(this->lookup("freeContactAngle"));
        }
    }

    // Check if correctPointNormals switch is set
    if (this->found("correctPointNormals"))
    {
        correctPointNormals_ = Switch(this->lookup("correctPointNormals"));
    }

    // Check if correctDisplacement switch is set
    if (this->found("correctDisplacement"))
    {
        correctDisplacement_ = Switch(this->lookup("correctDisplacement"));
    }

    // Check if correctCurvature switch is set
    if (this->found("correctCurvature"))
    {
        correctCurvature_ = Switch(this->lookup("correctCurvature"));
    }

    // Check if correctCurvature switch is set
    if (this->found("correctCurvature"))
    {
        smoothCurvature_ = Switch(this->lookup("smoothCurvature"));
    }

    // Check if curvExtrapOrder parameter is set
    if (this->found("curvExtrapOrder"))
    {
        curvExtrapOrder_ = Switch(this->lookup("curvExtrapOrder"));
    }

    // Check if fvcNGradUn switch is set
    if (this->found("fvcNGradUn"))
    {
        fvcNGradUn_ = Switch(this->lookup("fvcNGradUn"));
    }

    // Check if implicitCoupling switch is set
    if (this->found("implicitCoupling"))
    {
        implicitCoupling_ = Switch(this->lookup("implicitCoupling"));
    }

    // Check if interface deformation limit is set
    if (this->found("interfaceDeformationLimit"))
    {
        interfaceDeformationLimit_ =
            readScalar(this->lookup("interfaceDeformationLimit"));

        if (debug)
        {
            Info << "trackedSurface::updateProperties() : "
                << "Interface deformation limit: "
                << interfaceDeformationLimit_ << endl;
        }
    }
}


bool trackedSurface::updateMesh()
{
    if (!fixedInterface_ && totalDisplacementPtr_)
    {
        scalar minCellThickness =
            2*gMin(1.0/mesh().boundary()[aPatchID()].deltaCoeffs());

        scalar maxInterfaceDeformation =
            gMax(mag(totalDisplacement()))/minCellThickness;

        if (debug)
        {
            Info << "trackedSurface::updateMesh() : "
                << "Maximal relative interface deformation: "
                << maxInterfaceDeformation
                << endl;
        }

        // Move mesh only if interface deformation limit is exceeded
        if (maxInterfaceDeformation > interfaceDeformationLimit_)
        {
            if (debug)
            {
                Info << "trackedSurface::updateMesh() : "
                    << "Moving mesh."
                    << endl;
            }

            moveMeshPoints();

            moveUpdate();

            meshMoved_ = true;
        }
    }
    else
    {
        if (debug)
        {
            Info << "trackedSurface::updateMesh() : "
                << "Skipping mesh movement."
                << endl;
        }

        meshMoved_ = false;
    }

    if (!fixedInterface_)
    {
        resetControlPoints();

        adjustDisplacementDirections();
    }

    return meshMoved_;
}


void trackedSurface::updatePoints()
{
    scalarField& interfacePhi = phi().boundaryField()[aPatchID()];

    if (fixedInterface_)
    {
        interfacePhi = scalarField(interfacePhi.size(), 0.0);

        if (debug)
        {
            Info << "trackedSurface::updatePoints() : "
                << "Skipping point update."
                << endl;
        }
    }
    else
    {
        correctInterfacePhi(interfacePhi);

        for
        (
            int trackedSurfCorr=0;
            trackedSurfCorr<nTrackedSurfCorr_;
            trackedSurfCorr++
        )
        {
            movePoints(interfacePhi);
        }
    }
}


// TODO


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void trackedSurface::moveFixedPatches(const vectorField& displacement)
{
    // Take only displacement at fixed patches
    vectorField delta(aMesh().nPoints(), vector::zero);

    forAll (fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if (fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::moveFixedPatches(...): ")
                << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll (patchPoints, pointI)
        {
            delta[patchPoints[pointI]] = displacement[patchPoints[pointI]];
        }
    }

    moveMeshPoints(delta);

    moveUpdate();
}


// TODO
void trackedSurface::smoothing()
{
    if (smoothing_)
    {
    // Reset control points

        controlPoints() = aMesh().areaCentres().internalField();


    // Smoothing step 1

        scalarField deltaH = scalarField(controlPoints().size(), 0.0);

        vectorField displacement = pointDisplacement(deltaH);

        totalDisplacement() += displacement;

        if (total0Update_)
        {
            total0Displacement() += displacement;
        }

        vectorField newMeshPoints = calcNewMeshPoints(displacement);

        mesh().movePoints(newMeshPoints);


    // Smoothing step 2

        deltaH = calcDeltaH(-mesh().phi().boundaryField()[aPatchID()]
            *DB().deltaT().value());

        displacement = pointDisplacement(deltaH);

        totalDisplacement() += displacement;

        if (total0Update_)
        {
            total0Displacement() += displacement;
        }

        newMeshPoints = calcNewMeshPoints(displacement);

        mesh().movePoints(newMeshPoints);

    // Update all dependencies after point movement

        moveUpdate();
    }
}


// TODO
void trackedSurface::smoothMesh()
{
    const vectorField& oldPoints = aMesh().patch().localPoints();

    vectorField newPoints = oldPoints;

    const labelListList& pointEdges = aMesh().patch().pointEdges();

    const labelListList& pointFaces = aMesh().patch().pointFaces();

    const edgeList& edges = aMesh().patch().edges();

    const faceList& faces = aMesh().patch().localFaces();

    const labelList& boundaryPoints = aMesh().patch().boundaryPoints();

    Info << "Smoothing surface mesh" << endl;

    // Average edge length

    scalar avgEdgeLength = 0;

    forAll (edges, edgeI)
    {
        avgEdgeLength += edges[edgeI].mag(oldPoints);
    }
    avgEdgeLength /= edges.size();


    // Smooth boundary points

    forAll (aMesh().boundary(), patchI)
    {
        const labelList& pPointLabels =
            aMesh().boundary()[patchI].pointLabels();

        const labelListList& pPointEdges =
            aMesh().boundary()[patchI].pointEdges();

        const edgeList::subList pEdges =
            aMesh().boundary()[patchI].patchSlice(aMesh().edges());


        // Find fixed points
        boolList fixedPoints(pPointLabels.size(), false);

        forAll (fixedPoints, pointI)
        {
            if (pPointEdges[pointI].size() == 1)
            {
                fixedPoints[pointI] = true;
            }
        }


        // Perform smoothing
        scalarField residual(pPointLabels.size(), 0);
        label counter = 0;
        do
        {
            counter++;

            forAll (pPointLabels, pointI)
            {
                if (!fixedPoints[pointI])
                {
                    vector curNewPoint = vector::zero;

                    forAll (pPointEdges[pointI], eI)
                    {
                        label curEdgeIndex = pPointEdges[pointI][eI];

                        const edge& curEdge = pEdges[curEdgeIndex];

                        vector d =
                            newPoints
                            [
                                curEdge.otherVertex(pPointLabels[pointI])
                            ]
                          - newPoints[pPointLabels[pointI]];

                        curNewPoint += d;
                    }

                    curNewPoint /= pPointEdges[pointI].size();

                    curNewPoint += newPoints[pPointLabels[pointI]];


                    // Project new point to the interface

                    label nearestPointID = -1;
                    scalar minDist = GREAT;
                    forAll (pPointLabels, pI)
                    {
                        label curPoint = pPointLabels[pI];

                        scalar dist = mag(curNewPoint - oldPoints[curPoint]);

                        if (dist < minDist)
                        {
                            nearestPointID = pI;
                            minDist = dist;
                        }
                    }

                    bool foundProjection = false;
                    forAll (pPointEdges[nearestPointID], edgeI)
                    {
                        label edgeID = pPointEdges[nearestPointID][edgeI];

                        vector eTilda =
                            oldPoints
                            [
                                pEdges[edgeID]
                               .otherVertex(pPointLabels[nearestPointID])
                            ]
                          - oldPoints[pPointLabels[nearestPointID]];

                        scalar eMag = mag(eTilda);

                        eTilda /= mag(eTilda);

                        vector t =
                            eTilda
                           *(
                                eTilda
                               &(
                                    curNewPoint
                                  - oldPoints[pPointLabels[nearestPointID]]
                                )
                            );

                        if
                        (
                            ((t&eTilda) >= 0)
                         && ((t&eTilda) <= eMag)
                        )
                        {
                            curNewPoint =
                                oldPoints[pPointLabels[nearestPointID]] + t;
                            foundProjection = true;
                            break;
                        }
                    }

                    if (!foundProjection)
                    {
                        FatalErrorIn("trackedSurface::smoothMesh() : ")
                            << "Could not project patch point to surface"
                                << abort(FatalError);
                    }

                    residual[pointI] =
                        mag(curNewPoint - newPoints[pPointLabels[pointI]])
                       /(
                           mag(curNewPoint - oldPoints[pPointLabels[pointI]])
                         + SMALL
                        );

                    newPoints[pPointLabels[pointI]] = curNewPoint;
                }
            }
        }
        while(max(residual) > 1e-3);

        Info << "Patch: " << aMesh().boundary()[patchI].name()
            << ", max residual: " << max(residual)
            << ", num of iterations: " << counter << endl;
    }


    // Smooth internal points

    boolList fixedPoints(newPoints.size(), false);

    forAll (boundaryPoints, pointI)
    {
        fixedPoints[boundaryPoints[pointI]] = true;
    }

    scalarField residual(newPoints.size(), 0);
    label counter = 0;
    do
    {
        counter++;

        forAll (newPoints, pointI)
        {
            if (!fixedPoints[pointI])
            {
                vector curNewPoint = vector::zero;

                scalar sumW = 0;

                forAll (pointEdges[pointI], eI)
                {
                    label curEdgeIndex = pointEdges[pointI][eI];

                    const edge& curEdge = edges[curEdgeIndex];

                    vector d =
                        newPoints[curEdge.otherVertex(pointI)]
                      - newPoints[pointI];

//                     scalar w = 1.0;
                    scalar w = curEdge.mag(newPoints)/avgEdgeLength;

                    curNewPoint += w*d;

                    sumW += w;
                }

                curNewPoint /= sumW;

                curNewPoint += newPoints[pointI];


                // Project new point to the interface

                label nearestPointID = -1;
                scalar minDist = GREAT;
                forAll (oldPoints, pI)
                {
                    scalar dist = mag(curNewPoint - oldPoints[pI]);

                    if (dist < minDist)
                    {
                        nearestPointID = pI;
                        minDist = dist;
                    }
                }

                const vector& n = aMesh().pointAreaNormals()[nearestPointID];

                pointHit ph(curNewPoint);

                forAll (pointFaces[nearestPointID], faceI)
                {
                    label faceID = pointFaces[nearestPointID][faceI];

                    ph = faces[faceID].ray(curNewPoint, n, oldPoints);

                    if (ph.hit())
                    {
                        curNewPoint = ph.hitPoint();
                        break;
                    }
                }

                if (!ph.hit())
                {
                    Info << counter << ", " << pointI << endl;

                    FatalErrorIn("trackedSurface::smoothMesh() : ")
                        << "Could not project point to surface"
                            << abort(FatalError);
                }

                residual[pointI] =
                    mag(curNewPoint - newPoints[pointI])
                   /(
                        mag(curNewPoint - oldPoints[pointI])
                      + SMALL
                    );

                newPoints[pointI] = curNewPoint;
            }
        }
    }
    while(max(residual) > 1e-3);

    Info << "Internal points, max residual: " << max(residual)
            << ", num of iterations: " << counter << endl;

    // Mesh statistic

    scalar minEdge = GREAT;
    scalar avgEdge = 0;
    scalar maxEdge = SMALL;

    forAll (edges, edgeI)
    {
        scalar curEdgeLength = edges[edgeI].mag(newPoints);

        if (curEdgeLength < minEdge)
        {
            minEdge = curEdgeLength;
        }

        if (curEdgeLength > maxEdge)
        {
            maxEdge = curEdgeLength;
        }

        avgEdge += curEdgeLength;
    }
    avgEdge /= edges.size();

    Info << "Edge length, min: " << minEdge
        << ", max: " << maxEdge << ", avg: " << avgEdge << endl;


    vectorField displacement = newPoints - oldPoints;

    moveMeshPoints(displacement);

    moveUpdate();
}


// TODO
scalar trackedSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if (cleanInterface())
    {
        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/(cleanInterfaceSurfTension().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE =
            linearEdgeInterpolate(surfaceTension())().internalField()
          + SMALL;

        const scalarField& dE =aMesh().lPN();

        CoNum = gMax
        (
            DB().deltaT().value()/
            sqrt
            (
                rhoFluidA().value()*dE*dE*dE/
                2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


// TODO
void trackedSurface::smoothCurvature()
{
    areaScalarField& oldK =
        const_cast<areaScalarField&>
        (
// TEST: Sub-mesh
            curvature()
        );

    areaScalarField K
    (
        IOobject
        (
            "K",
            DB().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh(),
        dimless/dimLength,
        fixedGradientFaPatchScalarField::typeName
    );

    K.internalField() = oldK.internalField();

    forAll (K.boundaryField(), patchI)
    {
        if
        (
            K.boundaryField()[patchI].type()
         == fixedGradientFaPatchScalarField::typeName
        )
        {
            fixedGradientFaPatchScalarField& Kp =
                refCast<fixedGradientFaPatchScalarField>
                (
                    K.boundaryField()[patchI]
                );

            Kp.gradient() = 0;
        }
    }

    K.correctBoundaryConditions();


    label counter = 0;

    // Set gradient
    do
    {
        counter++;

        areaVectorField gradK = fac::grad(K);

        forAll (K.boundaryField(), patchI)
        {
            if
            (
                K.boundaryField()[patchI].type()
             == fixedGradientFaPatchScalarField::typeName
            )
            {
                fixedGradientFaPatchScalarField& Kp =
                    refCast<fixedGradientFaPatchScalarField>
                    (
                        K.boundaryField()[patchI]
                    );

                Kp.gradient() =
                (
                    aMesh().boundary()[patchI].edgeNormals()
                   &gradK.boundaryField()[patchI].patchInternalField()
                );
            }
        }

        K.correctBoundaryConditions();
    }
    while(counter<5);


    areaScalarField indicator =
        fac::div
        (
            fac::lnGrad(K)*aMesh().magLe()
          - (aMesh().Le()&fac::interpolate(fac::grad(K)))
        );

    scalar minIndicator = GREAT;
    label refFace = -1;

    forAll (indicator, faceI)
    {
        if (mag(indicator[faceI]) < minIndicator)
        {
            minIndicator = mag(indicator[faceI]);
            refFace = faceI;
        }
    }

    scalar gMinIndicator =
        returnReduce<scalar>(minIndicator, minOp<scalar>());

    bool procHasRef = false;
    if (mag(minIndicator - gMinIndicator) < SMALL)
    {
        procHasRef = true;
    }

    label procID = Pstream::nProcs();
    if (procHasRef)
    {
        procID = Pstream::myProcNo();
    }

    label minProcID =
        returnReduce<label>(procID, minOp<label>());

    if (procID != minProcID)
    {
        procHasRef = false;
    }

//     scalar refK = K[refFace];

    counter = 0;

    do
    {
        counter++;

        faScalarMatrix KEqn
        (
            fam::laplacian(K)
         == fac::div(aMesh().Le()&fac::interpolate(fac::grad(K)))
        );

//         KEqn.setReference(refFace, refK);
        if (K.needReference() && procHasRef)
        {
            KEqn.source()[refFace] +=
                KEqn.diag()[refFace]*K[refFace];

            KEqn.diag()[refFace] +=
                KEqn.diag()[refFace];
        }
        KEqn.solve();
    }
    while(counter<2);

    oldK = K;
}


// TODO
void trackedSurface::correctCurvature()
{
    // Correct curvature next to fixed patches

    areaScalarField& K =
        const_cast<areaScalarField&>
        (
// TEST: Sub-mesh
            curvature()
        );

    scalarField& KI = K.internalField();

    if (curvExtrapOrder_ == 0)
    {
        forAll (fixedTrackedSurfacePatches_, patchI)
        {
            label fixedPatchID =
                aMesh().boundary().findPatchID
                (
                    fixedTrackedSurfacePatches_[patchI]
                );

            if (fixedPatchID == -1)
            {
                FatalErrorIn("trackedSurface::trackedSurface(...)")
                    << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                        << " defined in the trackedSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            const labelListList& fFaces = aMesh().patch().faceFaces();

            forAll (eFaces, edgeI)
            {
                const label& curFace = eFaces[edgeI];
                const labelList& curFaceFaces = fFaces[curFace];

                scalar avrK = 0.0;
                label counter = 0;

                forAll (curFaceFaces, faceI)
                {
                    label index = findIndex(eFaces, curFaceFaces[faceI]);

                    if (index == -1)
                    {
                        avrK += K[curFaceFaces[faceI]];
                        counter++;
                    }
                }

                avrK /= counter;

                KI[curFace] = avrK;
            }
        }
    }
    else if (curvExtrapOrder_ == 1)
    {
        label counter = 0;
        do
        {
            counter++;

            K.correctBoundaryConditions();
            areaVectorField gradK = fac::grad(K);
            vectorField& gradKI = gradK.internalField();

            forAll (fixedTrackedSurfacePatches_, patchI)
            {
                label fixedPatchID =
                    aMesh().boundary().findPatchID
                    (
                        fixedTrackedSurfacePatches_[patchI]
                    );

                if (fixedPatchID == -1)
                {
                    FatalErrorIn("trackedSurface::trackedSurface(...)")
                        << "Wrong faPatch name in the fixedTrackedSurfacePatches "
                            << "list defined in the trackedSurfaceProperties "
                            << "dictionary"
                            << abort(FatalError);
                }

                const labelList& eFaces =
                    aMesh().boundary()[fixedPatchID].edgeFaces();

                const labelListList& fFaces = aMesh().patch().faceFaces();

                const vectorField& fCentres = aMesh().areaCentres();

                forAll (eFaces, edgeI)
                {
                    const label& curFace = eFaces[edgeI];
                    const labelList& curFaceFaces = fFaces[curFace];

                    scalar avrK = 0.0;
                    label counter = 0;

                    forAll (curFaceFaces, faceI)
                    {
                        label index = findIndex(eFaces, curFaceFaces[faceI]);

                        if (index == -1)
                        {
                            vector dr =
                                fCentres[curFace]
                              - fCentres[curFaceFaces[faceI]];

                            avrK += KI[curFaceFaces[faceI]]
                                + (dr&gradKI[curFaceFaces[faceI]]);
                            counter++;
                        }
                    }

                    avrK /= counter;

                    KI[curFace] = avrK;
                }
            }
        }
        while(counter<10);
    }
    else
    {
        FatalErrorIn("trackedSurface::trackedSurface(...)")
            << "Curvature extrapolation order is not set correctly"
                << abort(FatalError);
    }

    K.correctBoundaryConditions();
}


// TODO
void trackedSurface::correctPointNormals()
{
    // Correct normals for fixed patches points

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        );

    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const labelListList& fFaces =
        aMesh().patch().faceFaces();

    const faceList& faces =
        aMesh().patch().localFaces();

    const pointField& points =
        aMesh().patch().localPoints();


    // Wedge points
    forAll (aMesh().boundary(), patchI)
    {
        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPoints = wedgePatch.pointLabels();

            forAll (patchPoints, pointI)
            {
                label curPoint = patchPoints[pointI];

                labelHashSet faceSet;
                forAll (pFaces[curPoint], faceI)
                {
                    faceSet.insert(pFaces[curPoint][faceI]);
                }
                labelList curFaces = faceSet.toc();

                labelHashSet pointSet;

                pointSet.insert(curPoint);
                for(label i=0; i<curFaces.size(); i++)
                {
                    const labelList& facePoints = faces[curFaces[i]];
                    for(label j=0; j<facePoints.size(); j++)
                    {
                        if (!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();


                labelHashSet addPointsSet;
                forAll (curPoints, pointI)
                {
                    label index =
                        findIndex(patchPoints, curPoints[pointI]);

                    if (index != -1)
                    {
                        addPointsSet.insert(curPoints[pointI]);
                    }
                }
                addPointsSet.insert(curPoint);
                labelList curAddPoints = addPointsSet.toc();


                if (curPoints.size() + curAddPoints.size() >= 5)
                {
                    vectorField allPoints
                    (
                        curPoints.size()+curAddPoints.size()
                    );
                    scalarField W(curPoints.size()+curAddPoints.size(), 1.0);
                    label I = -1;
                    for(label i=0; i<curPoints.size(); i++)
                    {
                        I++;
                        allPoints[I] = points[curPoints[i]];
                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
                    }
                    for(label i=0; i<curAddPoints.size(); i++)
                    {
                        I++;
                        allPoints[I] =
                            transform
                            (
                                wedgePatch.faceT(),
                                points[curAddPoints[i]]
                            );
                        W[I] = 1.0/magSqr(allPoints[I] - points[curPoint]);
                    }

                    // Transforme points
                    vector origin = points[curPoint];
                    vector axis = N[curPoint]/mag(N[curPoint]);
                    vector dir = (allPoints[0] - points[curPoint]);
                    dir -= axis*(axis&dir);
                    dir /= mag(dir);
                    coordinateSystem cs("cs", origin, axis, dir);

                    forAll (allPoints, pI)
                    {
                        allPoints[pI] = cs.localPosition(allPoints[pI]);
                    }

                    scalarRectangularMatrix M
                    (
                        allPoints.size(),
                        5,
                        0.0
                    );

                    for(label i = 0; i < allPoints.size(); i++)
                    {
                        M[i][0] = sqr(allPoints[i].x());
                        M[i][1] = sqr(allPoints[i].y());
                        M[i][2] = allPoints[i].x()*allPoints[i].y();
                        M[i][3] = allPoints[i].x();
                        M[i][4] = allPoints[i].y();
                    }

                    scalarSquareMatrix MtM(5, 0.0);

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

                    scalarField MtR(5, 0);

                    for (label i=0; i<MtR.size(); i++)
                    {
                        for (label j=0; j<M.n(); j++)
                        {
                            MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                        }
                    }

                    scalarSquareMatrix::LUsolve(MtM, MtR);

                    vector curNormal = vector(MtR[3], MtR[4], -1);

                    curNormal = cs.globalVector(curNormal);

                    curNormal *= sign(curNormal&N[curPoint]);

                    N[curPoint] = curNormal;
                }
            }
        }
    }

    // Fixed boundary points

    forAll (fixedTrackedSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedTrackedSurfacePatches_[patchI]
            );

        if (fixedPatchID == -1)
        {
            FatalErrorIn("trackedSurface::trackedSurface(...)")
                << "Wrong faPatch name in the fixedTrackedSurfacePatches list"
                    << " defined in the trackedSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll (pLabels, pointI)
        {
            label curPoint = pLabels[pointI];

            labelHashSet faceSet;
            forAll (pFaces[curPoint], faceI)
            {
                faceSet.insert(pFaces[curPoint][faceI]);
            }

            labelList curFaces = faceSet.toc();

            forAll (curFaces, faceI)
            {
                const labelList& curFaceFaces =
                    fFaces[curFaces[faceI]];

                forAll (curFaceFaces, fI)
                {
                    label curFaceFace = curFaceFaces[fI];

                    label index = findIndex(eFaces, curFaceFace);

                    if ( (index==-1) && !faceSet.found(curFaceFace) )
                    {
                        faceSet.insert(curFaceFace);
                    }
                }
            }
            curFaces = faceSet.toc();

            labelHashSet pointSet;

            pointSet.insert(curPoint);
            for(label i=0; i<curFaces.size(); i++)
            {
                const labelList& fPoints = faces[curFaces[i]];
                for(label j=0; j<fPoints.size(); j++)
                {
                    if (!pointSet.found(fPoints[j]))
                    {
                        pointSet.insert(fPoints[j]);
                    }
                }
            }

            pointSet.erase(curPoint);

            labelList curPoints = pointSet.toc();

            // LS quadric fit
            vectorField allPoints(curPoints.size());
            scalarField W(curPoints.size(), 1.0);
            for(label i=0; i<curPoints.size(); i++)
            {
                allPoints[i] = points[curPoints[i]];
                W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
            }

            // Transforme points
            vector origin = points[curPoint];
            vector axis = N[curPoint]/mag(N[curPoint]);
            vector dir = (allPoints[0] - points[curPoint]);
            dir -= axis*(axis&dir);
            dir /= mag(dir);
            coordinateSystem cs("cs", origin, axis, dir);

            forAll (allPoints, pI)
            {
                allPoints[pI] = cs.localPosition(allPoints[pI]);
            }

            scalarRectangularMatrix M
            (
                allPoints.size(),
                5,
                0.0
            );

            for(label i = 0; i < allPoints.size(); i++)
            {
                M[i][0] = sqr(allPoints[i].x());
                M[i][1] = sqr(allPoints[i].y());
                M[i][2] = allPoints[i].x()*allPoints[i].y();
                M[i][3] = allPoints[i].x();
                M[i][4] = allPoints[i].y();
            }

            scalarSquareMatrix MtM(5, 0.0);

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

            scalarField MtR(5, 0);

            for (label i=0; i<MtR.size(); i++)
            {
                for (label j=0; j<M.n(); j++)
                {
                    MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                }
            }

            scalarSquareMatrix::LUsolve(MtM, MtR);

            vector curNormal = vector(MtR[3], MtR[4], -1);

            curNormal = cs.globalVector(curNormal);

            curNormal *= sign(curNormal&N[curPoint]);

            N[curPoint] = curNormal;
        }
    }

    // Correct wedge points
    forAll (aMesh().boundary(), patchI)
    {
        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPoints = wedgePatch.pointLabels();

            vector n =
                transform
                (
                    wedgePatch.edgeT(),
                    wedgePatch.centreNormal()
                );

            n /= mag(n);

            forAll (patchPoints, pointI)
            {
                N[patchPoints[pointI]]
                    -= n*(n&N[patchPoints[pointI]]);
            }
        }
    }


    // Boundary points correction
    forAll (aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().correctPatchPointNormals(patchI)
        && !aMesh().boundary()[patchI].coupled()
        )
        {
            if (aMesh().boundary()[patchI].ngbPolyPatchIndex() == -1)
            {
                FatalErrorIn
                    (
                        "void trackedSurface::correctPointNormals const"
                    )   << "Neighbour polyPatch index is not defined "
                        << "for faPatch " << aMesh().boundary()[patchI].name()
                        << abort(FatalError);
            }

            labelList patchPoints = aMesh().boundary()[patchI].pointLabels();

            vectorField n =
                aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

            forAll (patchPoints, pointI)
            {
                N[patchPoints[pointI]]
                    -= n[pointI]*(n[pointI]&N[patchPoints[pointI]]);
            }
        }
    }


    N /= mag(N);
}


// TODO: For steep contact angles, the code sometimes gets a fp exeption. Why?
// TODO: Is the definition of the contact angle right?
void trackedSurface::correctContactLinePointNormals()
{
    if (!freeContactAngle_)
    {
        // Correct normals for contact line points
        // according to specified contact angle

        vectorField& N =
            const_cast<vectorField&>
            (
                aMesh().pointAreaNormals()
            );

        if (contactAnglePtr_)
        {
            if (debug)
            {
                Info << "trackedSurface::correctContactLinePointNormals() : "
                    << "Correcting contact line normals."
                        << endl;
            }

            vectorField oldPoints(aMesh().nPoints(), vector::zero);

            const labelList& meshPoints = aMesh().patch().meshPoints();

            forAll (oldPoints, ptI)
            {
                oldPoints[ptI] =
                    mesh().oldPoints()[meshPoints[ptI]];
            }

#           include "createTangentField.H"

            forAll (aMesh().boundary(), patchI)
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

                        vectorField ngbN =
                            aMesh().boundary()[patchI].ngbPolyPatchPointNormals();

                        const labelList& patchPoints =
                            aMesh().boundary()[patchI].pointLabels();

                        vectorField pN = vectorField(N, patchPoints);

                        vectorField rotationAxis = (ngbN^pN);
                        rotationAxis /= mag(rotationAxis) + SMALL;


                        // Calc rotation axis using edge vectors

                        const edgeList& edges = aMesh().edges();

                        const labelListList& pointEdges =
                            aMesh().boundary()[patchI].pointEdges();

                        forAll (pointEdges, pointI)
                        {
                            vector rotAx = vector::zero;

                            vector rotAxisI = vector::zero;
                            scalar rotAngleI = 0;
                            scalar rotAngleWI = 0;
                            scalar rotAngleWsumI = 0;

                            forAll (pointEdges[pointI], edgeI)
                            {
                                label curEdge = pointEdges[pointI][edgeI];

                                label curGlobEdge =
                                    aMesh().boundary()[patchI].start() + curEdge;

                                vector e = edges[curGlobEdge].vec(oldPoints);

                                e *= (e&rotationAxis[pointI])
                                /mag(e&rotationAxis[pointI]);

                                e /= mag(e) + SMALL;

                                rotAx += e;

                                // Weight as inverse of edge lengths
                                rotAngleWI = 1.0/mag(e);

                                // Sum up weighted rotation angles and weights
                                rotAngleI += rotAngleWI * rotAngle[curEdge];
                                rotAngleWsumI += rotAngleWI;
                            }

                            if (pointEdges[pointI].size() == 1)
                            {
#                               include "addNgbProcessorEdgeTangent.H"
                            }

                            rotAxisI = rotAx/(mag(rotAx) + SMALL);

                            rotAngleI /= rotAngleWsumI;

                            vector oldNgbNI = ngbN[pointI];

                            // Rodrigues' rotation formula
                            ngbN[pointI] = oldNgbNI*cos(rotAngleI)
                            + rotAxisI*(rotAxisI & oldNgbNI)*(1 - cos(rotAngleI))
                            + (rotAxisI^oldNgbNI)*sin(rotAngleI);
                        }

                        forAll (patchPoints, pointI)
                        {
                            N[patchPoints[pointI]] -=
                                ngbN[pointI]*(ngbN[pointI]&N[patchPoints[pointI]]);

                            N[patchPoints[pointI]] /= mag(N[patchPoints[pointI]]);
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
