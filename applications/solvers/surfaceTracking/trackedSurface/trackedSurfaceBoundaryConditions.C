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

#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include "fixedGradientCorrectedFvPatchFields.H"
#include "fixedValueCorrectedFvPatchFields.H"

#include "zeroGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "mixedFvPatchFields.H"
#include "directionMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void trackedSurface::updateBoundaryConditions()
{
    if (!implicitCoupling_)
    {
        updateContactAngle();
        updateMuEff();
        updateSurfaceFlux();
        updateSurfactantConcentration();
        updateTemperature();
        updateConcentration();
        updateVelocity();
        updatePressure();
    }
}


void trackedSurface::updateContactAngle()
{
    // Correct contact angle acording to
    // current face normals

    if (freeContactAngle_)
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
                    // Calculate contact angle
                    scalarField& contactAngle =
                        contactAnglePtr_->boundaryField()[patchI];

                    const vectorField& nA =
                        aMesh().faceAreaNormals().internalField();

                    const vectorField ngbNA =
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals();

                    const unallocLabelList& edgeFaces =
                        aMesh().boundary()[patchI].edgeFaces();

                    forAll (edgeFaces, edgeI)
                    {
                        label faceI = edgeFaces[edgeI];

                        vector nAI = nA[faceI];
                        vector ngbNAI = ngbNA[edgeI];

                        contactAngle[edgeI] =
                            180.0/M_PI
                            * acos(-ngbNAI&nAI/mag(ngbNAI)/mag(nAI));
                    }
                }
            }
        }
    }
}


void trackedSurface::updateMuEff()
{
    if (turbulencePtr_)
    {
        const volScalarField nuEff = turbulence().nuEff();

        muEffFluidAval() =
            nuEff.boundaryField()[aPatchID()]
          * rho().boundaryField()[aPatchID()];

        if (twoFluids())
        {
            muEffFluidBval() =
                nuEff.boundaryField()[bPatchID()]
              * rho().boundaryField()[bPatchID()];
        }
    }
}


void trackedSurface::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void trackedSurface::updateSurfactantConcentration()
{
    if (!cleanInterface())
    {
        if (debug)
        {
            Info << "trackedSurface::updateSurfactantConcentration() : "
                << "Correct surfactant concentration."
                << endl;
        }

        // Crate and solve the surfactant transport equation
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


        if (surfactant().soluble())
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
                dimensionedScalar("Cb", dimMoles/dimVolume, 0),
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

        if (debug)
        {
            Info << "trackedSurface::updateSurfactantConcentration() : "
                << "Correct surface tension."
                << endl;
        }

        surfaceTension() =
            cleanInterfaceSurfTension()
          + surfactant().surfactR()
           *surfactant().surfactT()
           *surfactant().surfactSaturatedConc()
           *log(1.0 - surfactantConcentration()
           /surfactant().surfactSaturatedConc());

        if (neg(min(surfaceTension().internalField())))
        {
            FatalErrorIn("trackedSurface::correctSurfactantConcentration()")
                << "Surface tension is negative."
                    << abort(FatalError);
        }

        deleteDemandDrivenData(surfaceTensionForcePtr_);
    }
}


void trackedSurface::updateTemperature()
{
    if (TPtr_)
    {
        if (twoFluids())
        {
            // Update fixedValue boundary condition on patch B

            if
            (
                T().boundaryField()[bPatchID()].type()
             == fixedValueFvPatchScalarField::typeName
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
                FatalErrorIn("trackedSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for trackedSurfaceShadow patch is "
                    << T().boundaryField()[bPatchID()].type()
                    << ", instead of "
                    << fixedValueFvPatchScalarField::typeName
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
             == fixedGradientFvPatchScalarField::typeName
            )
            {
                fixedGradientFvPatchScalarField& aT =
                    refCast<fixedGradientFvPatchScalarField >
                    (
                        T().boundaryField()[aPatchID()]
                    );

                aT.gradient() = nGradT;
            }
            else
            {
                FatalErrorIn("trackedSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for trackedSurface patch is "
                    << T().boundaryField()[aPatchID()].type()
                    << ", instead of "
                    << fixedGradientFvPatchScalarField::typeName
                    << abort(FatalError);
            }
        }
        else
        {
            if
            (
                T().boundaryField()[aPatchID()].type()
             != zeroGradientFvPatchScalarField::typeName
             && T().boundaryField()[aPatchID()].type()
             != fixedGradientFvPatchScalarField::typeName
            )
            {
                FatalErrorIn("trackedSurface::updateTemperature()")
                    << "Bounary condition on " << T().name()
                    <<  " for trackedSurface patch is "
                    << T().boundaryField()[aPatchID()].type()
                    << ", instead of "
                    << zeroGradientFvPatchScalarField::typeName
                    << " or "
                    << fixedGradientFvPatchScalarField::typeName
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


void trackedSurface::updateConcentration()
{
    if (cPtr_)
    {
        if (twoFluids())
        {
            // Update fixedValue boundary condition on patch B

            if
            (
                c().boundaryField()[bPatchID()].type()
             == fixedValueFvPatchScalarField::typeName
            )
            {
                c().boundaryField()[bPatchID()] ==
                    interpolatorAB().faceInterpolate
                    (
                        c().boundaryField()[aPatchID()]
                    );
            }
            else
            {
                FatalErrorIn("trackedSurface::updateConcentration()")
                    << "Bounary condition on " << c().name()
                    <<  " for trackedSurfaceShadow patch is "
                    << c().boundaryField()[bPatchID()].type()
                    << ", instead of "
                    << fixedValueFvPatchScalarField::typeName
                    << abort(FatalError);
            }


            // Update fixedGradient boundary condition on patch A

            scalarField DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            scalarField cPB = interpolatorBA().faceInterpolate
            (
                c().boundaryField()[bPatchID()].patchInternalField()
            );

            const scalarField& cFs =
                c().boundaryField()[aPatchID()];

            scalarField nGradT = kFluidB().value()*(cPB - cFs)*DnA;
            nGradT /= kFluidA().value() + VSMALL;

            if
            (
                c().boundaryField()[aPatchID()].type()
             == fixedGradientFvPatchScalarField::typeName
            )
            {
                fixedGradientFvPatchScalarField& ac =
                    refCast<fixedGradientFvPatchScalarField >
                    (
                        c().boundaryField()[aPatchID()]
                    );

                ac.gradient() = nGradT;
            }
            else
            {
                FatalErrorIn("trackedSurface::updateConcentration()")
                    << "Bounary condition on " << c().name()
                    <<  " for trackedSurface patch is "
                    << c().boundaryField()[aPatchID()].type()
                    << ", instead of "
                    << fixedGradientFvPatchScalarField::typeName
                    << abort(FatalError);
            }
        }
        else
        {
            if
            (
                c().boundaryField()[aPatchID()].type()
             != zeroGradientFvPatchScalarField::typeName
             && c().boundaryField()[aPatchID()].type()
             != fixedGradientFvPatchScalarField::typeName
            )
            {
                FatalErrorIn("trackedSurface::updateConcentration()")
                    << "Bounary condition on " << c().name()
                    <<  " for trackedSurface patch is "
                    << c().boundaryField()[aPatchID()].type()
                    << ", instead of "
                    << zeroGradientFvPatchScalarField::typeName
                    << " or "
                    << fixedGradientFvPatchScalarField::typeName
                    << abort(FatalError);
            }
        }


        // Update surface tension

        concentration().internalField() =
            c().boundaryField()[aPatchID()];
        concentration().correctBoundaryConditions();

        const dictionary& concentrationDict =
            transport().subDict("concentration");

        dimensionedScalar ddcSigma(concentrationDict.lookup("ddcSigma"));
        dimensionedScalar cRefSigma(concentrationDict.lookup("cRefSigma"));

        surfaceTension() =
            cleanInterfaceSurfTension()
          + ddcSigma*(concentration() - cRefSigma);

        deleteDemandDrivenData(surfaceTensionForcePtr_);
    }
}


void trackedSurface::updateVelocity()
{
    const vectorField& nA = aMesh().faceAreaNormals().internalField();

    if (twoFluids())
    {
        vectorField nAf = mesh().boundary()[aPatchID()].nf();

        vectorField nBf = mesh().boundary()[bPatchID()].nf();

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
         == fixedGradientCorrectedFvPatchVectorField::typeName
        )
        {
            fixedGradientCorrectedFvPatchVectorField& aU =
                refCast<fixedGradientCorrectedFvPatchVectorField >
                (
                    U().boundaryField()[aPatchID()]
                );

            UPA += aU.corrVecGrad();
        }

        vectorField UtPA = UPA - nAf*(nAf & UPA);


        vectorField UPB = interpolatorBA().faceInterpolate
        (
            U().boundaryField()[bPatchID()].patchInternalField()
        );

        if
        (
            U().boundaryField()[bPatchID()].type()
         == fixedValueCorrectedFvPatchVectorField::typeName
        )
        {
            fixedValueCorrectedFvPatchVectorField& bU =
                refCast<fixedValueCorrectedFvPatchVectorField >
                (
                    U().boundaryField()[bPatchID()]
                );

            UPB += interpolatorBA().faceInterpolate(bU.corrVecGrad());
        }

        vectorField UtPB = UPB - nAf*(nAf & UPB);

        vectorField UtFs =
            muEffFluidAval()*DnA*UtPA
          + muEffFluidBval()*DnB*UtPB;

        vectorField UnFs =
            nA*phi().boundaryField()[aPatchID()]
           /mesh().boundary()[aPatchID()].magSf();

// // TEST: Use velocities to correct interface normal velocity
//         // Normal component
//         vectorField UnPA = nAf*(nAf & UPA);
//         vectorField UnPB = nAf*(nAf & UPB);
//
//         vectorField UnFs =
//             2*muEffFluidAval()*UnPA*DnA
//           + 2*muEffFluidBval()*UnPB*DnB;
//
//         UnFs /= 2*muEffFluidAval()*DnA
//           + 2*muEffFluidBval()*DnB + VSMALL;

        Us().internalField() += UnFs - nAf*(nAf&Us().internalField());
        Us().correctBoundaryConditions();

//*******************************************************************

        UtFs -= (muEffFluidAval() - muEffFluidBval())*
            (fac::grad(Us())().internalField()&nA);

        vectorField tangentialSurfaceTensionForce = ((I-nA*nA)&surfaceTensionForce());

        UtFs += tangentialSurfaceTensionForce;

        UtFs /= muEffFluidAval()*DnA + muEffFluidBval()*DnB + VSMALL;

        Us().internalField() = UnFs + UtFs;
        Us().correctBoundaryConditions();

        // Store old-time velocity field U()
        U().oldTime();

        U().boundaryField()[bPatchID()] ==
            interpolatorAB().faceInterpolate(UtFs)
          + nBf*fvc::meshPhi(rho(),U())().boundaryField()[bPatchID()]/
            mesh().boundary()[bPatchID()].magSf();

        if
        (
            p().boundaryField()[bPatchID()].type()
         == fixedGradientFvPatchScalarField::typeName
        )
        {
            fixedGradientFvPatchScalarField& pB =
                refCast<fixedGradientFvPatchScalarField >
                (
                    p().boundaryField()[bPatchID()]
                );

            pB.gradient() =
               - rhoFluidB().value()
                *(
                     nBf&fvc::ddt(U())().boundaryField()[bPatchID()]
                 );
        }


        // Update fixedGradient boundary condition on patch A

        updateNGradUn();

        vectorField nGradU =
            muEffFluidBval()*(UtPB - UtFs)*DnB // ZT, DnA
          + tangentialSurfaceTensionForce
          + muEffFluidAval()*nA*nGradUn()
          + (muEffFluidBval() - muEffFluidAval())
           *(fac::grad(Us())().internalField()&nA);

        nGradU /= muEffFluidAval() + VSMALL;


        if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientCorrectedFvPatchVectorField::typeName
        )
        {
            fixedGradientCorrectedFvPatchVectorField& aU =
                refCast<fixedGradientCorrectedFvPatchVectorField >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else if
        (
            U().boundaryField()[aPatchID()].type()
         == fixedGradientFvPatchVectorField::typeName
        )
        {
            fixedGradientFvPatchVectorField& aU =
                refCast<fixedGradientFvPatchVectorField >
                (
                    U().boundaryField()[aPatchID()]
                );

            aU.gradient() = nGradU;
        }
        else
        {
            FatalErrorIn("trackedSurface::updateVelocity()")
                << "Bounary condition on " << U().name()
                    <<  " for trackedSurface patch is "
                    << U().boundaryField()[aPatchID()].type()
                    << ", instead of "
                    << fixedGradientCorrectedFvPatchVectorField::typeName
                    << " or "
                    << fixedGradientFvPatchVectorField::typeName
                    << abort(FatalError);
        }
    }
    else
    {
        if (fixedInterface_)
        {
            vectorField tangentialSurfaceTensionForce =
                (I-sqr(nA)) & surfaceTensionForce();

            vectorField tnGradU =
                tangentialSurfaceTensionForce/(muEffFluidAval() + VSMALL);

            vectorField UPA =
                U().boundaryField()[aPatchID()].patchInternalField();

            const scalarField& DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            Us().internalField() = UPA + tnGradU/DnA;
            Us().internalField() -= sqr(nA) & Us().internalField();
            Us().correctBoundaryConditions();

            updateNGradUn();

            if
            (
                U().boundaryField()[aPatchID()].type()
             == directionMixedFvPatchVectorField::typeName
            )
            {
                directionMixedFvPatchVectorField& aU =
                    refCast<directionMixedFvPatchVectorField >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.refValue() = vector::zero;
                aU.refGrad() = tnGradU;
                aU.valueFraction() = sqr(nA);
            }
            else
            {
                FatalErrorIn("trackedSurface::updateVelocity()")
                    << "Bounary condition on " << U().name()
                        <<  " for trackedSurface patch is "
                        << U().boundaryField()[aPatchID()].type()
                        << ", instead of "
                        << directionMixedFvPatchVectorField::typeName
                        << abort(FatalError);
            }
        }
        else
        {
            vectorField UnFs =
                nA*phi_.boundaryField()[aPatchID()]
               /mesh().boundary()[aPatchID()].magSf();

// // TEST: Use velocity to correct interface normal velocity
//             vectorField UPA =
//                 U().boundaryField()[aPatchID()].patchInternalField();
//
//             if
//             (
//                 U().boundaryField()[aPatchID()].type()
//             == fixedGradientCorrectedFvPatchVectorField::typeName
//             )
//             {
//                 fixedGradientCorrectedFvPatchVectorField& aU =
//                     refCast<fixedGradientCorrectedFvPatchVectorField >
//                     (
//                         U().boundaryField()[aPatchID()]
//                     );
//
//                 UPA += aU.corrVecGrad();
//             }
//
//             vectorField UnFs = nA*(nA & UPA);

            // Correct normal component of surface velocity
            Us().internalField() += UnFs - nA*(nA&Us().internalField());
            Us().correctBoundaryConditions();

            vectorField tangentialSurfaceTensionForce = ((I-nA*nA)&surfaceTensionForce());

            vectorField tnGradU =
                tangentialSurfaceTensionForce/(muEffFluidAval() + VSMALL)
              - (fac::grad(Us())().internalField()&nA);

            vectorField UPA =
                U().boundaryField()[aPatchID()].patchInternalField();

            const scalarField& DnA = mesh().boundary()[aPatchID()].deltaCoeffs();

            vectorField UtFs = UPA + tnGradU/DnA;
            UtFs -= nA*(nA & UtFs);

            Us().internalField() = UtFs + UnFs;
            Us().correctBoundaryConditions();

            updateNGradUn();

            vectorField nGradU =
                tangentialSurfaceTensionForce/(muEffFluidAval() + VSMALL)
              + nA*nGradUn()
              - (fac::grad(Us())().internalField()&nA);

            if
            (
                U().boundaryField()[aPatchID()].type()
             == fixedGradientCorrectedFvPatchVectorField::typeName
            )
            {
                fixedGradientCorrectedFvPatchVectorField& aU =
                    refCast<fixedGradientCorrectedFvPatchVectorField >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.gradient() = nGradU;
            }
            else if
            (
                U().boundaryField()[aPatchID()].type()
             == fixedGradientFvPatchVectorField::typeName
            )
            {
                fixedGradientFvPatchVectorField& aU =
                    refCast<fixedGradientFvPatchVectorField >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.gradient() = nGradU;
            }
            else if
            (
                U().boundaryField()[aPatchID()].type()
             == mixedFvPatchVectorField::typeName
            )
            {
                mixedFvPatchVectorField& aU =
                    refCast<mixedFvPatchVectorField >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.refValue() = vector::zero;
                aU.refGrad() = nGradU;
                aU.valueFraction() = 0;
            }
            else if
            (
                U().boundaryField()[aPatchID()].type()
             == directionMixedFvPatchVectorField::typeName
            )
            {
                directionMixedFvPatchVectorField& aU =
                    refCast<directionMixedFvPatchVectorField >
                    (
                        U().boundaryField()[aPatchID()]
                    );

                aU.refValue() = vector::zero;
                aU.refGrad() = nGradU;
                aU.valueFraction() = symmTensor::zero;
            }
            else
            {
                FatalErrorIn("trackedSurface::updateVelocity()")
                    << "Bounary condition on " << U().name()
                        <<  " for trackedSurface patch is "
                        << U().boundaryField()[aPatchID()].type()
                        << ", instead of "
                        << fixedGradientCorrectedFvPatchVectorField::typeName
                        << ", "
                        << fixedGradientFvPatchVectorField::typeName
                        << " or "
                        << mixedFvPatchVectorField::typeName
                        << " or "
                        << directionMixedFvPatchVectorField::typeName
                        << abort(FatalError);
            }
        }
    }
}


void trackedSurface::updatePressure()
{
    const vectorField& nA = aMesh().faceAreaNormals().internalField();

    // Correct pressure boundary condition at the tracked surface

    if (twoFluids())
    {
        scalarField pA =
            interpolatorBA().faceInterpolate
            (
                p().boundaryField()[bPatchID()]
            );

        const scalarField& K = curvature().internalField();

        if (debug)
        {
            Info << "trackedSurface::updatePressure() : "
                << "Surface curvature: min = " << gMin(K)
                << ", max = " << gMax(K)
                << ", average = " << gAverage(K)
                << endl;
        }

        pA -= nA & surfaceTensionForce();

        pA += 2.0*(muEffFluidAval() - muEffFluidBval())*nGradUn();

        vector R0 = vector::zero;

        pA -= (rhoFluidA().value() - rhoFluidB().value())*
            (
                g_.value()
              & (
                    mesh().C().boundaryField()[aPatchID()]
                  - R0
                )
            );

        if (p0Ptr_)
        {
            pA += p0().boundaryField()[aPatchID()];
        }

        p().boundaryField()[aPatchID()] == pA;
    }
    else
    {
        if (fixedInterface_)
        {
            if
            (
               !(
                    p().boundaryField()[aPatchID()].type()
                 == fixedFluxPressureFvPatchScalarField::typeName
                 || p().boundaryField()[aPatchID()].type()
                 == zeroGradientFvPatchScalarField::typeName
                )
            )
            {
                FatalErrorIn("trackedSurface::updatePressure()")
                    << "Bounary condition on " << p().name()
                        <<  " for trackedSurface patch is "
                        << p().boundaryField()[aPatchID()].type()
                        << ", instead of "
                        << fixedFluxPressureFvPatchScalarField::typeName
                        << " or "
                        << zeroGradientFvPatchScalarField::typeName
                        << abort(FatalError);
            }
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

            const scalarField& K = curvature().internalField();

            if (debug)
            {
                Info << "trackedSurface::updatePressure() : "
                    << "Surface curvature: min = " << gMin(K)
                    << ", max = " << gMax(K)
                    << ", average = " << gAverage(K)
                    << endl;
            }

            pA -= nA & surfaceTensionForce();

            pA += 2.0*muEffFluidAval()*nGradUn();

            if (p0Ptr_)
            {
                pA += p0().boundaryField()[aPatchID()];
            }

            if
            (
               p().boundaryField()[aPatchID()].type()
            == fixedValueFvPatchScalarField::typeName
            )
            {
                fixedValueFvPatchScalarField& ap =
                    refCast<fixedValueFvPatchScalarField >
                    (
                        p().boundaryField()[aPatchID()]
                    );

                ap == pA;
            }
            else
            {
                FatalErrorIn("trackedSurface::updatePressure()")
                    << "Bounary condition on " << p().name()
                        <<  " for trackedSurface patch is "
                        << p().boundaryField()[aPatchID()].type()
                        << ", instead of "
                        << fixedValueFvPatchScalarField::typeName
                        << abort(FatalError);
            }
        }
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


void trackedSurface::updateNGradUn()
{
    if (fvcNGradUn_)
    {
        if (debug)
        {
            Info << "trackedSurface::updateNGradUn() : "
                << "Update normal derivative of normal velocity using fvc."
                << endl;
        }

        volVectorField phiU = fvc::reconstruct(phi());

        vectorField nA = mesh().boundary()[aPatchID()].nf();

        scalarField UnP =
            (nA&phiU.boundaryField()[aPatchID()].patchInternalField());

        scalarField UnFs =
            phi().boundaryField()[aPatchID()]
           /mesh().magSf().boundaryField()[aPatchID()];

        nGradUn() =
            (UnFs - UnP)*mesh().deltaCoeffs().boundaryField()[aPatchID()];

        bool secondOrderCorrection = true;

        if (secondOrderCorrection)
        {
            // Correct normal component of phiU
            // before gradient calculation
            forAll (phiU.boundaryField(), patchI)
            {
                vectorField n =
                    mesh().Sf().boundaryField()[patchI]
                   /mesh().magSf().boundaryField()[patchI];

                phiU.boundaryField()[patchI] +=
                    n
                   *(
                       (
                           phi().boundaryField()[patchI]
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
        if (debug)
        {
            Info << "trackedSurface::updateNGradUn() : "
                << "Update normal derivative of normal velocity using fac."
                << endl;
        }

        nGradUn() = -fac::div(Us())().internalField();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
