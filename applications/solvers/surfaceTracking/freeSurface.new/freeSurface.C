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

#include "freeSurface.H"

#include "volFields.H"
#include "transformField.H"

#include "emptyFaPatch.H"
#include "wedgeFaPatch.H"
#include "wallFvPatch.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroGradientCorrectedFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "primitivePatchInterpolation.H"

#include "coordinateSystem.H"
#include "scalarMatrices.H"
#include "fixedGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void freeSurface::clearOut()
{
    deleteDemandDrivenData(interpolatorABPtr_);
    deleteDemandDrivenData(interpolatorBAPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(totalDisplacementPtr_);
    deleteDemandDrivenData(aMeshPtr_);
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(fluidIndicatorPtr_);
    deleteDemandDrivenData(contactAnglePtr_);
    deleteDemandDrivenData(temperaturePtr_);
    deleteDemandDrivenData(surfaceTensionForcePtr_);
    deleteDemandDrivenData(nGradUnPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurface::freeSurface
(
    dynamicFvMesh& m,
    const volScalarField& rho,
    volVectorField& Ub,
    volScalarField& Pb,
    const surfaceScalarField& sfPhi,
    volScalarField* TbPtr
)
:
    IOdictionary
    (
        IOobject
        (
            "freeSurfaceProperties",
            Ub.mesh().time().constant(),
            Ub.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(m),
    rho_(rho),
    U_(Ub),
    p_(Pb),
    TPtr_(TbPtr),
    phi_(sfPhi),
    curTimeIndex_(Ub.mesh().time().timeIndex()),
    twoFluids_
    (
        this->lookup("twoFluids")
    ),
    normalMotionDir_
    (
        this->lookup("normalMotionDir")
    ),
    motionDir_(0, 0, 0),
    cleanInterface_
    (
        this->lookup("cleanInterface")
    ),
    aPatchID_(-1),
    bPatchID_(-1),
    muFluidA_
    (
        this->lookup("muFluidA")
    ),
    muFluidB_
    (
        this->lookup("muFluidB")
    ),
    rhoFluidA_
    (
        this->lookup("rhoFluidA")
    ),
    rhoFluidB_
    (
        this->lookup("rhoFluidB")
    ),
    kFluidA_("kFluidA", dimThermalConductivity, 0.0),
    kFluidB_("kFluidB", dimThermalConductivity, 0.0),
    CpFluidA_("CpFluidA", dimSpecificHeatCapacity, 0.0),
    CpFluidB_("CpFluidB", dimSpecificHeatCapacity, 0.0),
    g_(this->lookup("g")),
    cleanInterfaceSurfTension_
    (
        this->lookup("surfaceTension")
    ),
    fixedFreeSurfacePatches_
    (
        this->lookup("fixedFreeSurfacePatches")
    ),
    pointNormalsCorrectionPatches_
    (
        this->lookup("pointNormalsCorrectionPatches")
    ),
    nFreeSurfCorr_
    (
        readInt(this->lookup("nFreeSurfaceCorrectors"))
    ),
    smoothing_(false),
    correctPointNormals_(false),
    correctDisplacement_(false),
    correctCurvature_(false),
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
    aMeshPtr_(NULL),
    UsPtr_(NULL),
    phisPtr_(NULL),
    surfactConcPtr_(NULL),
    surfaceTensionPtr_(NULL),
    surfactantPtr_(NULL),
    fluidIndicatorPtr_(NULL),
    contactAnglePtr_(NULL),
    temperaturePtr_(NULL),
    surfaceTensionForcePtr_(NULL),
    nGradUnPtr_(NULL)
{
    //Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = vector(this->lookup("motionDir"));
        motionDir_ /= mag(motionDir_) + SMALL;
    }

    // Set point normal correction patches
    boolList& correction = aMesh().correctPatchPointNormals();

    forAll(pointNormalsCorrectionPatches_, patchI)
    {
        word patchName = pointNormalsCorrectionPatches_[patchI];

        label patchID = aMesh().boundary().findPatchID(patchName);

        if(patchID == -1)
        {
            FatalErrorIn
            (
                "freeSurface::freeSurface(...)"
            )   << "Patch name for point normals correction does not exist"
                << abort(FatalError);
        }

        correction[patchID] = true;
    }

    // Detect the free surface patch
    forAll (mesh().boundary(), patchI)
    {
        if(mesh().boundary()[patchI].name() == "freeSurface")
        {
            aPatchID_ = patchI;

            Info<< "Found free surface patch. ID: " << aPatchID_
                << endl;
        }
    }

    if(aPatchID() == -1)
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Free surface patch not defined.  Please make sure that "
                << " the free surface patches is named as freeSurface"
                << abort(FatalError);
    }


    // Detect the free surface shadow patch
    if (twoFluids())
    {
        forAll (mesh().boundary(), patchI)
        {
            if(mesh().boundary()[patchI].name() == "freeSurfaceShadow")
            {
                bPatchID_ = patchI;

                Info<< "Found free surface shadow patch. ID: "
                    << bPatchID_ << endl;
            }
        }

        if(bPatchID() == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Free surface shadow patch not defined. "
                    << "Please make sure that the free surface shadow patch "
                    << "is named as freeSurfaceShadow."
                    << abort(FatalError);
        }
    }

    // Mark free surface boundary points
    // which belonge to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }

    // Mark fixed free surface boundary points
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }

    // Mark free-surface boundary point
    // at the axis of 2-D axisymmetic cases
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            forAll(wedgePatch.axisPoints(), pI)
            {
                motionPointsMask()[wedgePatch.axisPoints()[pI]] = 0;
            }

//             if(wedgePatch.axisPoint() > -1)
//             {
//                 motionPointsMask()[wedgePatch.axisPoint()] = 0;

//                 Info << "Axis point: "
//                     << wedgePatch.axisPoint()
//                     << " vector: "
//                     << aMesh().points()[wedgePatch.axisPoint()] << endl;
//             }
        }
    }

    // Read free-surface points total displacement if present
    readTotalDisplacement();

    // Read control points positions if present
    controlPoints();

    // Check if smoothing switch is set
    if (this->found("smoothing"))
    {
        smoothing_ = Switch(this->lookup("smoothing"));
    }

    // Check if contactAngle is defined
    IOobject contactAngleHeader
    (
        "contactAngle",
        DB().timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (contactAngleHeader.headerOk())
    {
        Pout << "Reading contact angle" << endl;

        contactAnglePtr_ =
            new edgeScalarField
            (
                IOobject
                (
                    "contactAngle",
                    DB().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh()
            );
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

    // Check if curvExtrapOrder parameter is set
    if (this->found("curvExtrapOrder"))
    {
        curvExtrapOrder_ = Switch(this->lookup("curvExtrapOrder"));
    }

    // Check Marangoni effect
    if (TPtr_ && !cleanInterface())
    {
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Marangoni effect due to both "
                << "surfactant concentration gradient "
                << "and temperature gradient is not implemented"
                << abort(FatalError);
    }

    if (TPtr_)
    {
        // Read free-surface temperature field
        makeTemperature();

        // Read properties
        kFluidA_ = dimensionedScalar(this->lookup("kFluidA"));
        kFluidB_ = dimensionedScalar(this->lookup("kFluidB"));
        CpFluidA_ = dimensionedScalar(this->lookup("CpFluidA"));
        CpFluidB_ = dimensionedScalar(this->lookup("CpFluidB"));
    }

    // Clear geometry
    aMesh().movePoints();

    // Contact angle correction
    correctContactLinePointNormals();

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

        Info << "Interface deformation limit: "
            << interfaceDeformationLimit_ << endl;
    }
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

freeSurface::~freeSurface()
{
    clearOut();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void freeSurface::updateDisplacementDirections()
{
    if(normalMotionDir())
    {
        // Update point displacement correction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correcte point displacement direction
        // at the "centerline" symmetryPlane which represents the axis
        // of an axisymmetric case
        forAll(aMesh().boundary(), patchI)
        {
            if(aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
            {
                const wedgeFaPatch& wedgePatch =
                    refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

                vector axis = wedgePatch.axis();

                label centerLinePatchID =
                    aMesh().boundary().findPatchID("centerline");

                if(centerLinePatchID != -1)
                {
                    const labelList& pointLabels =
                        aMesh().boundary()[centerLinePatchID].pointLabels();

                    forAll(pointLabels, pointI)
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
                    Info << "Warning: centerline polyPatch does not exist. "
                        << "Free surface points displacement directions "
                        << "will not be corrected at the axis (centerline)"
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


bool freeSurface::predictPoints()
{
    // Smooth interface

    if (smoothing_)
    {
        smoothing();
    }

    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
//         moveMeshPoints();
    }

    return true;
}


bool freeSurface::correctPoints()
{
    for
    (
        int freeSurfCorr=0;
        freeSurfCorr<nFreeSurfCorr_;
        freeSurfCorr++
    )
    {
        movePoints(phi_.boundaryField()[aPatchID()]);
//         moveMeshPoints();
    }

    return true;
}


bool freeSurface::movePoints(const scalarField& interfacePhi)
{
    pointField newMeshPoints = mesh().allPoints();

    scalarField sweptVolCorr =
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

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
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    pointField displacement = pointDisplacement(deltaH);
//     displacement *= 0;

//     Pout << gMax(mag(displacement)) << endl;

    if (correctDisplacement_)
    {
//         Pout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
        correctPointDisplacement(sweptVolCorr, displacement);
    }

    // Move only free-surface points

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
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


    // Update total displacement field

//     if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
//     {
//         FatalErrorIn("freeSurface::movePoints()")
//             << "Total displacement of free surface points "
//                 << "from previous time step is not absorbed by the mesh."
//                 << abort(FatalError);
//     }
//     else if (curTimeIndex_ < DB().timeIndex())
//     {
//         totalDisplacement() = displacement;

//         curTimeIndex_ = DB().timeIndex();
//     }
//     else
//     {
//         totalDisplacement() += displacement;
//     }

    totalDisplacement() += displacement;

    twoDPointCorrector twoDPointCorr(mesh());

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints();

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//         correctCurvature();
        smoothCurvature();
    }

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
            ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& pU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            pU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
            ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& pP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

            pP.movePatchSubMesh();
        }
    }

    return true;
}


bool freeSurface::moveMeshPointsForOldFreeSurfDisplacement()
{
    if(totalDisplacementPtr_)
    {
        scalar minCellThickness =
            2*gMin(1.0/mesh().boundary()[aPatchID()].deltaCoeffs());

        scalar maxInterfaceDeformation =
            gMax(mag(totalDisplacement()))/minCellThickness;

        Info << "Maximal relative interface deformation: "
            << maxInterfaceDeformation << endl;


        // Move whole mesh only if interface deformation limit is exceeded
        if (maxInterfaceDeformation > interfaceDeformationLimit_)
        {
            Info << "Moving whole mesh" << endl;

            pointField newPoints = mesh().allPoints();

            const labelList& meshPointsA =
                mesh().boundaryMesh()[aPatchID()].meshPoints();

            forAll (totalDisplacement(), pointI)
            {
                newPoints[meshPointsA[pointI]] -= totalDisplacement()[pointI];
            }


            // Check mesh motion solver type

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
                        totalDisplacement()/DB().deltaT().value()
                    );

                if(twoFluids_)
                {
                    const labelList& meshPointsB =
                        mesh().boundaryMesh()[bPatchID()].meshPoints();

                    pointField totDisplacementB =
                        interpolatorAB().pointInterpolate
                        (
                            totalDisplacement()
                        );

                    forAll (totDisplacementB, pointI)
                    {
                        newPoints[meshPointsB[pointI]] -=
                            totDisplacementB[pointI];
                    }

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
                            totDisplacementB/DB().deltaT().value()
                        );
                }
            }
            else if (fvMotionSolver)
            {
                pointVectorField& motionU =
                    const_cast<pointVectorField&>
                    (
                        mesh().objectRegistry::
                        lookupObject<pointVectorField>
                        (
                            "pointMotionU"
                        )
                    );

                fixedValuePointPatchVectorField& motionUaPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[aPatchID()]
                    );

                motionUaPatch ==
                    totalDisplacement()/DB().deltaT().value();

                if(twoFluids_)
                {
                    const labelList& meshPointsB =
                        mesh().boundaryMesh()[bPatchID()].meshPoints();

                    pointField totDisplacementB =
                        interpolatorAB().pointInterpolate
                        (
                            totalDisplacement()
                        );

                    forAll (totDisplacementB, pointI)
                    {
                        newPoints[meshPointsB[pointI]] -=
                            totDisplacementB[pointI];
                    }

                    fixedValuePointPatchVectorField& motionUbPatch =
                        refCast<fixedValuePointPatchVectorField>
                        (
                            motionU.boundaryField()[bPatchID()]
                        );

                    motionUbPatch ==
                        totDisplacementB/DB().deltaT().value();
                }
            }

            twoDPointCorrector twoDPointCorr(mesh());

            twoDPointCorr.correctPoints(newPoints);

            mesh().movePoints(newPoints);

            deleteDemandDrivenData(totalDisplacementPtr_);

            mesh().update();

            aMesh().movePoints();

            correctContactLinePointNormals();

            if (correctPointNormals_)
            {
                correctPointNormals();
            }

            if (correctCurvature_)
            {
//             correctCurvature();
                smoothCurvature();
            }


            // Move correctedFvPatchField fvSubMeshes

            forAll(U().boundaryField(), patchI)
            {
                if
                (
                    (
                        U().boundaryField()[patchI].type()
                     == fixedGradientCorrectedFvPatchField<vector>::typeName
                    )
                ||
                    (
                        U().boundaryField()[patchI].type()
                     == fixedValueCorrectedFvPatchField<vector>::typeName
                    )
                ||
                    (
                        U().boundaryField()[patchI].type()
                     == zeroGradientCorrectedFvPatchField<vector>::typeName
                    )
                )
                {
                    correctedFvPatchField<vector>& aU =
                        refCast<correctedFvPatchField<vector> >
                        (
                            U().boundaryField()[patchI]
                        );

                    aU.movePatchSubMesh();
                }
            }

            forAll(p().boundaryField(), patchI)
            {
                if
                (
                    (
                        p().boundaryField()[patchI].type()
                     == fixedGradientCorrectedFvPatchField<scalar>::typeName
                    )
                ||
                    (
                        p().boundaryField()[patchI].type()
                     == fixedValueCorrectedFvPatchField<scalar>::typeName
                    )
                ||
                    (
                        p().boundaryField()[patchI].type()
                     == zeroGradientCorrectedFvPatchField<scalar>::typeName
                    )
                )
                {
                    correctedFvPatchField<scalar>& aP =
                        refCast<correctedFvPatchField<scalar> >
                        (
                            p().boundaryField()[patchI]
                        );

                    aP.movePatchSubMesh();
                }
            }
        }
    }

    return true;
}


bool freeSurface::smoothMesh()
{
    const vectorField& oldPoints = aMesh().patch().localPoints();

    vectorField newPoints = oldPoints;

    const labelListList& pointEdges = aMesh().patch().pointEdges();

    const labelListList& pointFaces = aMesh().patch().pointFaces();

    const edgeList& edges = aMesh().patch().edges();

    const faceList& faces = aMesh().patch().localFaces();

    const labelList& boundaryPoints = aMesh().patch().boundaryPoints();

    Info << "Smoothing free-surface mesh" << endl;

    // Average edge length

    scalar avgEdgeLength = 0;

    forAll(edges, edgeI)
    {
        avgEdgeLength += edges[edgeI].mag(oldPoints);
    }
    avgEdgeLength /= edges.size();


    // Smooth boundary points

    forAll(aMesh().boundary(), patchI)
    {
        const labelList& pPointLabels =
            aMesh().boundary()[patchI].pointLabels();

        const labelListList& pPointEdges =
            aMesh().boundary()[patchI].pointEdges();

        const edgeList::subList pEdges =
            aMesh().boundary()[patchI].patchSlice(aMesh().edges());


        // Find fixed points
        boolList fixedPoints(pPointLabels.size(), false);

        forAll(fixedPoints, pointI)
        {
            if (pPointEdges[pointI].size() == 1)
            {
                fixedPoints[pointI] = true;
            }
        }


        // Performe smoothing
        scalarField residual(pPointLabels.size(), 0);
        label counter = 0;
        do
        {
            counter++;

            forAll(pPointLabels, pointI)
            {
                if (!fixedPoints[pointI])
                {
                    vector curNewPoint = vector::zero;

                    forAll(pPointEdges[pointI], eI)
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
                    forAll(pPointLabels, pI)
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
                    forAll(pPointEdges[nearestPointID], edgeI)
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
                        FatalErrorIn("freeSurface::smoothMesh()")
                            << "Could not project patch point to free-surface"
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

    forAll(boundaryPoints, pointI)
    {
        fixedPoints[boundaryPoints[pointI]] = true;
    }

    scalarField residual(newPoints.size(), 0);
    label counter = 0;
    do
    {
        counter++;

        forAll(newPoints, pointI)
        {
            if (!fixedPoints[pointI])
            {
                vector curNewPoint = vector::zero;

                scalar sumW = 0;

                forAll(pointEdges[pointI], eI)
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
                forAll(oldPoints, pI)
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

                forAll(pointFaces[nearestPointID], faceI)
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

                    FatalErrorIn("freeSurface::smoothMesh()")
                        << "Could not project point to free-surface"
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

    forAll(edges, edgeI)
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

    //-- Set mesh motion boundary conditions

    // Check mesh motion solver type

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

        if(twoFluids_)
        {
            pointField displacementB =
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
    }
    else if (fvMotionSolver)
    {
        pointVectorField& motionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& motionUaPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        motionUaPatch ==
            displacement/DB().deltaT().value();

        if(twoFluids_)
        {
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );

            fixedValuePointPatchVectorField& motionUbPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            motionUbPatch ==
                displacementB/DB().deltaT().value();
        }
    }

    mesh().update();

    aMesh().movePoints();

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
        smoothCurvature();
    }

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& aU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            aU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& aP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

                aP.movePatchSubMesh();
        }
    }

    return true;
}


bool freeSurface::moveMeshPoints(const scalarField& interfacePhi)
{
    scalarField sweptVolCorr =
        interfacePhi
      - fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];

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
        sweptVolCorr *= (2.0/3.0)*DB().deltaT().value();
    }
    else
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Unsupported temporal differencing scheme : "
                << ddtScheme
                << abort(FatalError);
    }

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    scalarField deltaH =
        sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    pointField displacement = pointDisplacement(deltaH);

    //-- Set mesh motion boundary conditions

    // Check mesh motion solver type

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

        if(twoFluids_)
        {
            pointField displacementB =
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
    }
    else if (fvMotionSolver)
    {
        pointVectorField& motionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& motionUaPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        motionUaPatch ==
            displacement/DB().deltaT().value();

        if(twoFluids_)
        {
            pointField displacementB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );

            fixedValuePointPatchVectorField& motionUbPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            motionUbPatch ==
                displacementB/DB().deltaT().value();
        }
    }

    mesh().update();

    aMesh().movePoints();

    correctContactLinePointNormals();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//         correctCurvature();
        smoothCurvature();
    }

    // Move correctedFvPatchField fvSubMeshes

    forAll(U().boundaryField(), patchI)
    {
        if
        (
            (
                U().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<vector>::typeName
            )
         ||
            (
                U().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<vector>::typeName
            )
        )
        {
            correctedFvPatchField<vector>& aU =
                refCast<correctedFvPatchField<vector> >
                (
                    U().boundaryField()[patchI]
                );

            aU.movePatchSubMesh();
        }
    }

    forAll(p().boundaryField(), patchI)
    {
        if
        (
            (
                p().boundaryField()[patchI].type()
             == fixedGradientCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == fixedValueCorrectedFvPatchField<scalar>::typeName
            )
         ||
            (
                p().boundaryField()[patchI].type()
             == zeroGradientCorrectedFvPatchField<scalar>::typeName
            )
        )
        {
            correctedFvPatchField<scalar>& aP =
                refCast<correctedFvPatchField<scalar> >
                (
                    p().boundaryField()[patchI]
                );

                aP.movePatchSubMesh();
        }
    }

    return true;
}


void freeSurface::updateBoundaryConditions()
{
    if (!implicitCoupling_)
    {
        updateTemperature();
        updateVelocity();
        updateSurfactantConcentration();
        updatePressure();
    }
}


void freeSurface::updateVelocity()
{
    if(twoFluids())
    {
        vectorField nA = mesh().boundary()[aPatchID()].nf();

        vectorField nB = mesh().boundary()[bPatchID()].nf();

        scalarField DnB = interpolatorBA().faceInterpolate
        (
            mesh().boundary()[bPatchID()].deltaCoeffs()
        );

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();


        vectorField UPA =
            U().boundaryField()[aPatchID()].patchInternalField();

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            UPA += aU.corrVecGrad();
        }

        vectorField UtPA = UPA - nA*(nA & UPA);


        vectorField UPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedValueCorrectedFvPatchField<vector>& bU =
                refCast<fixedValueCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[bPatchID()]
                );

            UPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        vectorField UtPB = UPB - nA*(nA & UPB);

        vectorField UtFs =
            muFluidA().value()*DnA*UtPA
          + muFluidB().value()*DnB*UtPB;

        // Normal component
        vectorField UnPA = nA*(nA & UPA);
        vectorField UnPB = nA*(nA & UPB);

        vectorField UnFs =
            2*muFluidA().value()*UnPA*DnA
          + 2*muFluidB().value()*UnPB*DnB;

        UnFs /= 2*muFluidA().value()*DnA
          + 2*muFluidB().value()*DnB + VSMALL;

//         vectorField UnFs =
//             nA*phi_.boundaryField()[aPatchID()]
//            /mesh().boundary()[aPatchID()].magSf();

        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        correctUsBoundaryConditions();

//*******************************************************************

        UtFs -= (muFluidA().value() - muFluidB().value())*
            (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();


        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(MarangoniStress())
        {
            Info << "MarangoniStress------------------------" << endl;

            const vectorField& totSurfTensionForce = surfaceTensionForce();

            tangentialSurfaceTensionForce = ((I-nA*nA)&totSurfTensionForce);

//             tangentialSurfaceTensionForce =
//                 surfaceTensionGrad()().internalField();
        }
        else
        {
//             vectorField surfaceTensionForce =
//                 cleanInterfaceSurfTension().value()
//                *fac::edgeIntegrate
//                 (
//                     aMesh().Le()*aMesh().edgeLengthCorrection()
//                 )().internalField();

//             tangentialSurfaceTensionForce =
//                 surfaceTensionForce
//               - cleanInterfaceSurfTension().value()
//                *aMesh().faceCurvatures().internalField()*nA;
        }

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muFluidA().value()*DnA + muFluidB().value()*DnB + VSMALL;

        Us().internalField() = UnFs + UtFs;
        correctUsBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] ==
            interpolatorAB().faceInterpolate(UtFs)
          + nB*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchField<scalar>::typeName
        )
        {
            fixedGradientFvPatchField<scalar>& pB =
                refCast<fixedGradientFvPatchField<scalar> >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() =
               - rhoFluidB().value()
                *(
                     nB&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }


        // Update fixedGradient boundary condition on patch A

        updateNGradUn();

        vectorField nGradU =
            muFluidB().value()*(UtPB - UtFs)*DnB // ZT, DnA
          + tangentialSurfaceTensionForce
          + muFluidA().value()*nA*nGradUn()
          + (muFluidB().value() - muFluidA().value())
           *(fac::grad(Us())().internalField()&nA);

        nGradU /= muFluidA().value() + VSMALL;


        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for freeSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }
    else
    {
        const vectorField& nA = aMesh().faceAreaNormals().internalField();

        vectorField UnFs =
            nA*phi_.boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

        // Correct normal component of free-surface velocity
        Us().internalField() += UnFs - nA*(nA&Us().internalField());
        Us().correctBoundaryConditions();

        vectorField tangentialSurfaceTensionForce(nA.size(), vector::zero);

        if(MarangoniStress())
        {
            const vectorField& totSurfTensionForce = surfaceTensionForce();

            tangentialSurfaceTensionForce =
                ((I-nA*nA)&totSurfTensionForce);

//             tangentialSurfaceTensionForce =
//                 surfaceTensionGrad()().internalField();
        }
        else
        {
//             vectorField surfaceTensionForce =
//                 cleanInterfaceSurfTension().value()
//                *fac::edgeIntegrate
//                 (
//                     aMesh().Le()*aMesh().edgeLengthCorrection()
//                 )().internalField();

//             tangentialSurfaceTensionForce =
//                 surfaceTensionForce
//               - cleanInterfaceSurfTension().value()
//                *aMesh().faceCurvatures().internalField()*nA;

//             if (muFluidA().value() < SMALL)
//             {
//                 tangentialSurfaceTensionForce = vector::zero;
//             }
        }

        vectorField tnGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          - (fac::grad(Us())&aMesh().faceAreaNormals())().internalField();

        vectorField UtPA =
            U().boundaryField()[aPatchID()].patchInternalField();
        UtPA -= nA*(nA & UtPA);

        scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

        vectorField UtFs = UtPA + tnGradU/DnA;

        Us().internalField() = UtFs + UnFs;
        Us().correctBoundaryConditions();

        updateNGradUn();

//         scalarField divUs = -nGradUn();

        vectorField nGradU =
            tangentialSurfaceTensionForce/(muFluidA().value() + VSMALL)
          + nA*nGradUn()
          - (fac::grad(Us())().internalField()&nA);

        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchField<vector>::typeName
        )
        {
            fixedGradientCorrectedFvPatchField<vector>& aU =
                refCast<fixedGradientCorrectedFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchField<vector>::typeName
        )
        {
            fixedGradientFvPatchField<vector>& aU =
                refCast<fixedGradientFvPatchField<vector> >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("freeSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for freeSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientCorrectedFvPatchField<vector>::typeName
                    << " or "
                    << fixedGradientFvPatchField<vector>::typeName
                    << abort(FatalError);
        }
    }
}


void freeSurface::updatePressure()
{
    // Correct pressure boundary condition at the free-surface

    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if(twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        if(!MarangoniStress())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            if (correctCurvature_)
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;

                pA -= surfTensionK;
            }
            else
            {
                const vectorField& nA =
                    aMesh().faceAreaNormals().internalField();

                const vectorField& totSurfTensionForce = surfaceTensionForce();

                scalarField surfTensionK = (nA&totSurfTensionForce);

                pA -= surfTensionK;
            }
        }


//         fixedGradientFvPatchField<vector>& aU =
//             refCast<fixedGradientFvPatchField<vector> >
//             (
//                 U().boundaryField()[aPatchID()]
//             );

//         pA += 2.0*(muFluidA().value() - muFluidB().value())
//            *(nA&aU.gradient());

        pA += 2.0*(muFluidA().value() - muFluidB().value())*nGradUn();



//         pA -= 2.0*(muFluidA().value() - muFluidB().value())
//             *fac::div(Us())().internalField();

        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        vector R0 = vector::zero;

        scalarField pA =
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K)
            << endl;

        if(!MarangoniStress())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            if (correctCurvature_)
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;

                pA -= surfTensionK;
            }
            else
            {
                const vectorField& nA =
                    aMesh().faceAreaNormals().internalField();
                const vectorField& totSurfTensionForce = surfaceTensionForce();

                scalarField surfTensionK = (nA&totSurfTensionForce);

                pA -= surfTensionK;
            }
        }

//         scalarField divUs = -nGradUn();

        pA += 2.0*muFluidA().value()*nGradUn();
//         pA -= 2.0*muFluidA().value()*fac::div(Us())().internalField();

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

    vector R0 = vector::zero;

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::updatePressure(const scalarField& p0)
{
    // Correct pressure boundary condition at the free-surface

    vectorField nA = mesh().boundary()[aPatchID()].nf();

    if(twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

        if(!MarangoniStress())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            if (correctCurvature_)
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;

                pA -= surfTensionK;
            }
            else
            {
                const vectorField& nA =
                    aMesh().faceAreaNormals().internalField();

                const vectorField& totSurfTensionForce = surfaceTensionForce();

                scalarField surfTensionK = (nA&totSurfTensionForce);

                pA -= surfTensionK;
            }
        }

        pA += 2.0*(muFluidA().value() - muFluidB().value())*nGradUn();
//         pA += 2.0*(muFluidA().value() - muFluidB().value())
//             *fac::div(Us())().internalField();

        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        vector R0 = vector::zero;

        scalarField pA = p0
          - rhoFluidA().value()*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K)
            << endl;

        if(!MarangoniStress())
        {
//             pA -= cleanInterfaceSurfTension().value()*(K - gAverage(K));
            pA -= cleanInterfaceSurfTension().value()*K;
        }
        else
        {
            if (correctCurvature_)
            {
                scalarField surfTensionK =
                    surfaceTension().internalField()*K;

                pA -= surfTensionK;
            }
            else
            {
                const vectorField& nA =
                    aMesh().faceAreaNormals().internalField();
                const vectorField& totSurfTensionForce = surfaceTensionForce();

                scalarField surfTensionK = (nA&totSurfTensionForce);

                pA -= surfTensionK;
            }
        }

//         scalarField divUs = -nGradUn();

        pA += 2.0*muFluidA().value()*nGradUn();
//         pA += 2.0*muFluidA().value()*fac::div(Us())().internalField();

        p().boundaryField()[aPatchID()] == pA;
    }


    // Set modified pressure at patches with fixed apsolute
    // pressure

    vector R0 = vector::zero;

    for (int patchI=0; patchI < p().boundaryField().size(); patchI++)
    {
        if
        (
            p().boundaryField()[patchI].type()
         == fixedValueFvPatchScalarField::typeName
        )
        {
            if (patchI != aPatchID())
            {
                p().boundaryField()[patchI] ==
                  - rho().boundaryField()[patchI]
                   *(g_.value()&(mesh().C().boundaryField()[patchI] - R0));
            }
        }
    }
}


void freeSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void freeSurface::updateSurfactantConcentration()
{
    if(!cleanInterface())
    {
        Info << "Correct surfactant concentration" << endl << flush;

        updateSurfaceFlux();

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().surfactDiffusion(),
                surfactantConcentration()
            )
        );


        if(surfactant().soluble())
        {
            const scalarField& C =
                mesh().boundary()[aPatchID()]
               .lookupPatchField<volScalarField, scalar>("C");

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    DB().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensioned<scalar>("Cb", dimMoles/dimVolume, 0),
                zeroGradientFaPatchScalarField::typeName
            );

            Cb.internalField() = C;
            Cb.correctBoundaryConditions();

            CsEqn +=
                fam::Sp
                (
                    surfactant().surfactAdsorptionCoeff()*Cb
                  + surfactant().surfactAdsorptionCoeff()
                   *surfactant().surfactDesorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().surfactAdsorptionCoeff()
               *Cb*surfactant().surfactSaturatedConc();
        }

        CsEqn.solve();

        Info << "Correct surface tension" << endl;

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());

        if(neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn
            (
                "void freeSurface::correctSurfactantConcentration()"
            )
                << "Surface tension is negative"
                    << abort(FatalError);
        }
    }
}


void freeSurface::updateTemperature()
{
    if (TPtr_)
    {
        if(twoFluids())
        {
            // Update fixedValue boundary condition on patch B

            if
            (
                T().boundaryField()[bPatchID()].type()
             == fixedValueFvPatchField<scalar>::typeName
            )
            {
                T().boundaryField()[bPatchID()] ==
                    interpolatorAB().faceInterpolate
                    (
                        T().boundaryField()[aPatchID()]
                    );
            }
            else
            {
                FatalErrorIn("freeSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for freeSurfaceShadow patch is "
                    << T().boundaryField()[bPatchID()].type()
                    << ", instead "
                    << fixedValueFvPatchField<scalar>::typeName
                    << abort(FatalError);
            }


            // Update fixedGradient boundary condition on patch A

            scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            scalarField TPB = interpolatorBA().faceInterpolate
            (
                T().boundaryField()[bPatchID()].patchInternalField()
            );

            const scalarField& TFs =
                T().boundaryField()[aPatchID()];

            scalarField nGradT = kFluidB().value()*(TPB - TFs)*DnA;
            nGradT /= kFluidA().value() + VSMALL;

            if
            (
                T().boundaryField()[aPatchID()].type()
             == fixedGradientFvPatchField<scalar>::typeName
            )
            {
                fixedGradientFvPatchField<scalar>& aT =
                    refCast<fixedGradientFvPatchField<scalar> >
                    (
                        T().boundaryField()[aPatchID()]
                    );

                aT.gradient() = nGradT;
            }
            else
            {
                FatalErrorIn("freeSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for freeSurface patch is "
                    << T().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << fixedGradientFvPatchField<scalar>::typeName
                    << abort(FatalError);
            }
        }
        else
        {
            if
            (
                T().boundaryField()[aPatchID()].type()
             != zeroGradientFvPatchField<scalar>::typeName
            )
            {
                FatalErrorIn("freeSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for freeSurface patch is "
                    << T().boundaryField()[aPatchID()].type()
                    << ", instead "
                    << zeroGradientFvPatchField<scalar>::typeName
                    << abort(FatalError);
            }
        }


        // Update surface tension

        temperature().internalField() =
            T().boundaryField()[aPatchID()];
        temperature().correctBoundaryConditions();

        // Correct surface tension
        dimensionedScalar thermalCoeff
        (
            this->lookup("thermalCoeff")
        );

        dimensionedScalar refTemperature
        (
            this->lookup("refTemperature")
        );

        surfaceTension() =
            cleanInterfaceSurfTension()
          + thermalCoeff*(temperature() - refTemperature);

        deleteDemandDrivenData(surfaceTensionForcePtr_);
    }
}


void freeSurface::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryField()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if(ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U().boundaryField()[ngbPolyPatchID].type()
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type()
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N =
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


vector freeSurface::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[aPatchID()];

    vectorField pressureForces = S*P*n;

    return gSum(pressureForces);
}


vector freeSurface::totalViscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU =
        U().boundaryField()[aPatchID()].snGrad();

    vectorField viscousForces =
      - muFluidA().value()*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        );

    return gSum(viscousForces);
}


vector freeSurface::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), vector::zero);

    if(cleanInterface())
    {
        surfTensionForces =
            S*cleanInterfaceSurfTension().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces *= surfaceTension().internalField()*K;
    }

    return gSum(surfTensionForces);
}


void freeSurface::initializeControlPointsPosition()
{
    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();


    pointField newPoints = points + displacement;

    scalarField sweptVol(faces.size(), 0.0);

    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), vector::zero);

    forAll (faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].normal(newPoints);
    }

    forAll(deltaH, faceI)
    {
        deltaH[faceI] = sweptVol[faceI]/
            ((faceArea[faceI] & facesDisplacementDir()[faceI]) + SMALL);

        if ((faceArea[faceI] & facesDisplacementDir()[faceI]) < SMALL)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Something is probably wrong with the specified motion direction"
                    << abort(FatalError);

        }
    }

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    displacement = pointDisplacement(deltaH);
}


void freeSurface::smoothing()
{
    controlPoints() = aMesh().areaCentres().internalField();

    scalarField deltaH = scalarField(aMesh().nFaces(), 0.0);

    pointField displacement = pointDisplacement(deltaH);

    vectorField newMeshPoints = mesh().points();

    const labelList& meshPointsA =
        mesh().boundaryMesh()[aPatchID()].meshPoints();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorrector twoDPointCorr(mesh());
    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints();

    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    deltaH =
       -mesh().phi().boundaryField()[aPatchID()]*DB().deltaT().value()
       /(Sf*(Nf & facesDisplacementDir()));

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH[eFaces[edgeI]] *= 2.0;
        }
    }

    displacement = pointDisplacement(deltaH);

    newMeshPoints = mesh().points();

    forAll (displacement, pointI)
    {
        newMeshPoints[meshPointsA[pointI]] += displacement[pointI];
    }

    if(twoFluids_)
    {
        const labelList& meshPointsB =
            mesh().boundaryMesh()[bPatchID_].meshPoints();

        pointField displacementB =
            interpolatorAB().pointInterpolate
            (
                displacement
            );

        forAll (displacementB, pointI)
        {
            newMeshPoints[meshPointsB[pointI]] += displacementB[pointI];
        }
    }

    // Update total displacement field

    if(totalDisplacementPtr_ && (curTimeIndex_ < DB().timeIndex()))
    {
        FatalErrorIn("freeSurface::movePoints()")
            << "Total displacement of free surface points "
                << "from previous time step is not absorbed by the mesh."
                << abort(FatalError);
    }
    else if (curTimeIndex_ < DB().timeIndex())
    {
        totalDisplacement() = displacement;

        curTimeIndex_ = DB().timeIndex();
    }
    else
    {
        totalDisplacement() += displacement;
    }

    twoDPointCorr.correctPoints(newMeshPoints);

    mesh().movePoints(newMeshPoints);

    aMesh().movePoints();

    if (correctPointNormals_)
    {
        correctPointNormals();
    }

    if (correctCurvature_)
    {
//         correctCurvature();
        smoothCurvature();
    }
}


scalar freeSurface::maxCourantNumber()
{
    scalar CoNum = 0;

    if(cleanInterface())
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


bool freeSurface::MarangoniStress() const
{
    return (!cleanInterface() || TPtr_);
}


void freeSurface::updateProperties()
{
    muFluidA_ = dimensionedScalar(this->lookup("muFluidA"));

    muFluidB_ = dimensionedScalar(this->lookup("muFluidB"));

    rhoFluidA_ = dimensionedScalar(this->lookup("rhoFluidA"));

    rhoFluidB_ = dimensionedScalar(this->lookup("rhoFluidB"));

    g_ = dimensionedVector(this->lookup("g"));

    cleanInterfaceSurfTension_ =
        dimensionedScalar(this->lookup("surfaceTension"));

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

    // Check if curvExtrapOrder parameter is set
    if (this->found("curvExtrapOrder"))
    {
        curvExtrapOrder_ = Switch(this->lookup("curvExtrapOrder"));
    }
}


void freeSurface::writeVTK() const
{
    aMesh().patch().writeVTK
    (
        DB().timePath()/"freeSurface",
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void freeSurface::writeVTKControlPoints()
{
    // Write patch and points into VTK
    fileName name(DB().timePath()/"freeSurfaceControlPoints");
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << controlPoints().size() << " float" << nl;

    forAll(controlPoints(), pointI)
    {
        mps << controlPoints()[pointI].x() << ' '
            << controlPoints()[pointI].y() << ' '
            << controlPoints()[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << controlPoints().size() << ' '
        << controlPoints().size()*2 << nl;

    forAll(controlPoints(), pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


void freeSurface::correctCurvature()
{
    // Correct curvature next to fixed patches

    areaScalarField& K =
        const_cast<areaScalarField&>
        (
            aMesh().faceCurvatures()
        );

    scalarField& KI = K.internalField();

    if (curvExtrapOrder_ == 0)
    {
        forAll(fixedFreeSurfacePatches_, patchI)
        {
            label fixedPatchID =
                aMesh().boundary().findPatchID
                (
                    fixedFreeSurfacePatches_[patchI]
                );

            if(fixedPatchID == -1)
            {
                FatalErrorIn("freeSurface::freeSurface(...)")
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                        << " defined in the freeSurfaceProperties dictionary"
                        << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            const labelListList& fFaces = aMesh().patch().faceFaces();

            forAll(eFaces, edgeI)
            {
                const label& curFace = eFaces[edgeI];
                const labelList& curFaceFaces = fFaces[curFace];

                scalar avrK = 0.0;
                label counter = 0;

                forAll(curFaceFaces, faceI)
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

            forAll(fixedFreeSurfacePatches_, patchI)
            {
                label fixedPatchID =
                    aMesh().boundary().findPatchID
                    (
                        fixedFreeSurfacePatches_[patchI]
                    );

                if(fixedPatchID == -1)
                {
                    FatalErrorIn("freeSurface::freeSurface(...)")
                        << "Wrong faPatch name in the fixedFreeSurfacePatches "
                            << "list defined in the freeSurfaceProperties "
                            << "dictionary"
                            << abort(FatalError);
                }

                const labelList& eFaces =
                    aMesh().boundary()[fixedPatchID].edgeFaces();

                const labelListList& fFaces = aMesh().patch().faceFaces();

                const vectorField& fCentres = aMesh().areaCentres();

                forAll(eFaces, edgeI)
                {
                    const label& curFace = eFaces[edgeI];
                    const labelList& curFaceFaces = fFaces[curFace];

                    scalar avrK = 0.0;
                    label counter = 0;

                    forAll(curFaceFaces, faceI)
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
        FatalErrorIn("freeSurface::freeSurface(...)")
            << "Curvature extrapolation order is not set correctly"
                << abort(FatalError);
    }

    K.correctBoundaryConditions();
}


void freeSurface::correctPointNormals()
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
    forAll(aMesh().boundary(), patchI)
    {
        if (aMesh().boundary()[patchI].type() == wedgeFaPatch::typeName)
        {
            const wedgeFaPatch& wedgePatch =
                refCast<const wedgeFaPatch>(aMesh().boundary()[patchI]);

            const labelList& patchPoints = wedgePatch.pointLabels();

            forAll(patchPoints, pointI)
            {
                label curPoint = patchPoints[pointI];

                labelHashSet faceSet;
                forAll(pFaces[curPoint], faceI)
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
                        if(!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();


                labelHashSet addPointsSet;
                forAll(curPoints, pointI)
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

                    forAll(allPoints, pI)
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

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(pLabels, pointI)
        {
            label curPoint = pLabels[pointI];

            labelHashSet faceSet;
            forAll(pFaces[curPoint], faceI)
            {
                faceSet.insert(pFaces[curPoint][faceI]);
            }

            labelList curFaces = faceSet.toc();

            forAll(curFaces, faceI)
            {
                const labelList& curFaceFaces =
                    fFaces[curFaces[faceI]];

                forAll(curFaceFaces, fI)
                {
                    label curFaceFace = curFaceFaces[fI];

                    label index = findIndex(eFaces, curFaceFace);

                    if( (index==-1) && !faceSet.found(curFaceFace) )
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
                    if(!pointSet.found(fPoints[j]))
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

            forAll(allPoints, pI)
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

    // Correcte wedge points
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
                        "void freeSurface::correctPointNormals const"
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


void freeSurface::correctContactLinePointNormals()
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
        Pout << "Correcting contact line normals" << endl;

        vectorField oldPoints(aMesh().nPoints(), vector::zero);

        const labelList& meshPoints = aMesh().patch().meshPoints();

        forAll(oldPoints, ptI)
        {
            oldPoints[ptI] =
                mesh().oldPoints()[meshPoints[ptI]];
        }

#       include "createTangentField.H"

        forAll(aMesh().boundary(), patchI)
        {
            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    scalar rotAngle =
                        average
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        );

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

                        forAll(pointEdges[pointI], eI)
                        {
                            label curEdge =
                                aMesh().boundary()[patchI].start()
                              + pointEdges[pointI][eI];

                            vector e = edges[curEdge].vec(oldPoints);

                            e *= (e&rotationAxis[pointI])
                               /mag(e&rotationAxis[pointI]);

                            e /= mag(e) + SMALL;

                            rotAx += e;
                        }

                        if (pointEdges[pointI].size() == 1)
                        {
#                           include "addNgbProcessorEdgeTangent.H"
                        }

                        rotationAxis[pointI] = rotAx/(mag(rotAx) + SMALL);
                    }

                    // Rodrigues' rotation formula
                    ngbN = ngbN*cos(rotAngle)
                      + rotationAxis*(rotationAxis & ngbN)*(1 - cos(rotAngle))
                      + (rotationAxis^ngbN)*sin(rotAngle);

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


void freeSurface::correctPointDisplacement
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

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::freeSurface(...)")
                 << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        labelHashSet pointSet;

        forAll(eFaces, edgeI)
        {
            label curFace = eFaces[edgeI];

            const labelList& curPoints = faces[curFace];

            forAll(curPoints, pointI)
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

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;

            forAll(pFaces[curPoint], faceI)
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

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll(curPointFaces, faceI)
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


void freeSurface::smoothCurvature()
{
    areaScalarField& oldK =
        const_cast<areaScalarField&>
        (
            aMesh().faceCurvatures()
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

    forAll(K.boundaryField(), patchI)
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

        forAll(K.boundaryField(), patchI)
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

    forAll(indicator, faceI)
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


void freeSurface::updateNGradUn()
{
    if (fvcNGradUn_)
    {
        Info << "Update normal derivative of normal velocity using fvc"
            << endl;

        volVectorField phiU = fvc::reconstruct(phi_);

        vectorField nA = mesh().boundary()[aPatchID()].nf();

        scalarField UnP =
            (nA&phiU.boundaryField()[aPatchID()].patchInternalField());

        scalarField UnFs =
            phi_.boundaryField()[aPatchID()]
           /mesh().magSf().boundaryField()[aPatchID()];

        nGradUn() =
            (UnFs - UnP)*mesh().deltaCoeffs().boundaryField()[aPatchID()];

        bool secondOrderCorrection = true;

        if (secondOrderCorrection)
        {
            // Correct normal component of phiU
            // befor gradient calculation
            forAll(phiU.boundaryField(), patchI)
            {
                vectorField n =
                    mesh().Sf().boundaryField()[patchI]
                   /mesh().magSf().boundaryField()[patchI];

                phiU.boundaryField()[patchI] +=
                    n
                   *(
                       (
                           phi_.boundaryField()[patchI]
                          /mesh().magSf().boundaryField()[patchI]
                       )
                     - (n&phiU.boundaryField()[patchI])
                    );
            }

            // Calc gradient
            tensorField gradPhiUp =
                fvc::grad(phiU)().boundaryField()[aPatchID()]
               .patchInternalField();

            nGradUn() = 2*nGradUn() - (nA&(gradPhiUp&nA));
        }
    }
    else
    {
        Info << "Update normal derivative of normal velocity using fac"
            << endl;

        nGradUn() = -fac::div(Us())().internalField();

//         Us().correctBoundaryConditions();

//         areaVectorField UsTmp = Us();

//         vector avgUs = gAverage(UsTmp.internalField());
//         Pout << avgUs << endl;
//         UsTmp = dimensionedVector("avgUs", Us().dimensions(), avgUs);
//         UsTmp.correctBoundaryConditions();

//         edgeVectorField eUs = fac::interpolate(Us());
//         eUs = dimensionedVector("avgUs", Us().dimensions(), vector(1,0,0));

//         Pout << aMesh().weights().boundaryField() << endl;




//         nGradUn() =
//            -fac::edgeIntegrate
//             (
//                 aMesh().Le() & fac::interpolate(Us())
//             )
//           + (fac::edgeIntegrate(aMesh().Le()) & Us());

//         nGradUn() *= 2;




//           + aMesh().faceCurvatures()*(aMesh().faceAreaNormals()&Us());
    }

//     nGradUn() *= 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
