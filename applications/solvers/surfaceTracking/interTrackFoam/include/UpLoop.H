

    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UpLoop.H : "
                << "Commencing PIMPLE U-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            uniformDimensionedVectorField& g = storage.g();
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
// TODO: Make trackedSurface use real rho and mu fields!
            volScalarField& rho = storage.rho();
            twoPhaseMixture& transport = storage.transport();
            incompressible::turbulenceModel& turbulence = storage.turbulence();
            trackedSurface& interface = storage.interface();

            if (settings.relToUinf)
            {
                // Update free-stream force
// TODO: Check force calculation!
//                 storage.fInf().value() =
//                     interface.totalPressureForce();
//                 storage.fInf().value() =
//                     interface.totalViscousForce();
                storage.fInf().value() =
                    interface.totalPressureForce()
                  + interface.totalViscousForce();
//                 storage.fInf().value() =
//                     interface.totalForce();
//                 storage.fInf().value() =
//                   interface.totalNormalForce();
//                 storage.fInf().value() =
//                   interface.totalTangentialForce();
//                 storage.fInf().value() =
//                   interface.totalNormalForce()
//                 - interface.totalNormalSurfaceTensionForce();

                dimensionedVector Finf
                (
                    "Finf",
                    dimForce,
                    vector::zero
                );

                forAll(mesh.boundary(), patchI)
                {
                    if
                    (
                        mesh.boundary()[patchI].type()
                     != wedgeFvPatch::typeName
                    )
                    {
                        const vectorField& Sf = mesh.Sf().boundaryField()[patchI];

                        const scalarField& P = p.boundaryField()[patchI];

                        symmTensorField tau =
                            (
                                rho*turbulence.nu()
                            * dev(twoSymm(fvc::grad(U)))
                            )().boundaryField()[patchI];
//                           - rho.boundaryField()[patchI]
//                           * turbulence.devReff()().boundaryField()[patchI];

                        Finf.value() += gSum(Sf & (P*I - tau));
                    }
                }

                dimensionedVector FinfVol
                (
                    "Finf",
                    fvc::domainIntegrate
                    (
                        fvc::div
                        (
                            p*I
                          - rho*turbulence.nu()*dev(twoSymm(fvc::grad(U)))
                        )
                    )
                );

                // Constraint free-stream force
                forAll(mesh.boundary(), patchI)
                {
//                     if
//                     (
//                         mesh.boundary()[patchI].type()
//                      == wedgeFvPatch::typeName
//                     )
//                     {
//                         const wedgePolyPatch& wp =
//                             refCast<const wedgePolyPatch>
//                             (
//                                 mesh.boundary()[patchI].patch()
//                             );
//
//                         vector ez = wp.axis();
//                         vector ephi = wp.centreNormal();
//                         vector er = ez ^ ephi;
//
//                         storage.fInf().value() -=
//                             sqr(er) & storage.fInf().value();
//                     }

// TODO: Hard coded!!!
                    storage.fInf().value() -=
                        sqr(vector(1,0,0)) & storage.fInf().value();
                    storage.fInf().value() -=
                        sqr(vector(0,0,1)) & storage.fInf().value();

                    Finf.value() -=
                        sqr(vector(1,0,0)) & Finf.value();
                    Finf.value() -=
                        sqr(vector(0,0,1)) & Finf.value();

                    FinfVol.value() -=
                        sqr(vector(1,0,0)) & FinfVol.value();
                    FinfVol.value() -=
                        sqr(vector(0,0,1)) & FinfVol.value();
                }

                // Store old free-stream velocity
                dimensionedVector Uinf(storage.Uinf());
                if (control.firstIter())
                {
                    storage.UinfOld() = storage.Uinf();
                }

                // Init U
                if (max(mag(U)).value() == 0)
                {
                    U = Uinf;
                }

                // Calculate acceleration from du = F/m * dt
                dimensionedVector aInf
                (
                    "aInf",
                    storage.fInf() / storage.mInf()

                );

//                 // Calculate acceleration from du = F/m * dt
//                 dimensionedVector aInf
//                 (
//                     "aInf",
//                     FinfVol / fvc::domainIntegrate(rho)
//
//                 );

                // Calculate change of Uinf
                dimensionedVector dUinf
                (
                    "dUinf",
                   - aInf * runTime.deltaT()
                );

                // Update free-stream velocity with relaxation
                storage.Uinf() =
                    storage.UinfOld() + storage.dUinfRelax() * dUinf;

                Info<< "Free-stream velocity adjustment:"
                    << " fInf = " << storage.fInf().value()
                    << " , Uinf = " << storage.Uinf().value()
                    << " , UinfOld = " << storage.UinfOld().value()
                    << " , deltaT = " << runTime.deltaT().value()
                    << " , aInf = " << aInf.value()
                    << " , dUinf = " << dUinf.value()
                    << " , max(mag(U)) = " << max(mag(U)).value()
                    << " , Finf = " << Finf.value()
                    << " , FinfVol = " << FinfVol.value()
                    << endl;

                // Adjust boundary conditions for velocity
                {
                    forAll(U.boundaryField(), patchI)
                    {
                        if
                        (
                            U.boundaryField()[patchI].type()
                        == fixedValueFvPatchVectorField::typeName
                        && patchI != interface.aPatchID()
                        )
                        {
                            U.boundaryField()[patchI] ==
                                storage.Uinf().value();
                        }
                    }

                    U.correctBoundaryConditions();
                }

                // Adjust c
                storage.c() += storage.gradcInf() & mesh.C();

                // Adjust boundary conditions for c
                {
                    forAll(storage.c().boundaryField(), patchI)
                    {
                        if
                        (
                            storage.c().boundaryField()[patchI].type()
                         == fixedGradientFvPatchScalarField::typeName
                         && patchI == interface.aPatchID()
                        )
                        {
                            fixedGradientFvPatchScalarField& pc =
                                refCast<fixedGradientFvPatchScalarField >
                                (
                                    storage.c().boundaryField()[interface.aPatchID()]
                                );

                            pc.gradient() = 0.0;
                        }
                        else if
                        (
                            storage.c().boundaryField()[patchI].type()
                         != fixedGradientFvPatchScalarField::typeName
                         && patchI == interface.aPatchID()
                        )
                        {
                            FatalErrorIn(Control::typeName + " | UpLoop.H")
                                << "Bounary condition on " << storage.c().name()
                                    <<  " for trackedSurface patch is "
                                    << storage.c().boundaryField()[interface.aPatchID()].type()
                                    << ", instead of "
                                    << fixedGradientFvPatchScalarField::typeName
                                    << abort(FatalError);
                        }
                    }

                    storage.c().correctBoundaryConditions();
                }
            }

            interface.updateBoundaryConditions(); // abs phi

            if (settings.relToUinf)
            {
                // Adjust c
                storage.c() -= storage.gradcInf() & mesh.C();

                // Adjust boundary conditions for c
                {
                    forAll(storage.c().boundaryField(), patchI)
                    {
                        if
                        (
                            storage.c().boundaryField()[patchI].type()
                         == fixedGradientFvPatchScalarField::typeName
                         && patchI == interface.aPatchID()
                        )
                        {
                            fixedGradientFvPatchScalarField& pc =
                                refCast<fixedGradientFvPatchScalarField >
                                (
                                    storage.c().boundaryField()[interface.aPatchID()]
                                );

                            const vectorField n =
                                mesh.Sf().boundaryField()[patchI]
                              / mesh.magSf().boundaryField()[patchI];

                            pc.gradient() = - n & storage.gradcInf().value();
                        }
                        else if
                        (
                            storage.c().boundaryField()[patchI].type()
                         != fixedGradientFvPatchScalarField::typeName
                         && patchI == interface.aPatchID()
                        )
                        {
                            FatalErrorIn(Control::typeName + " | UpLoop.H")
                                << "Bounary condition on " << storage.c().name()
                                    <<  " for trackedSurface patch is "
                                    << storage.c().boundaryField()[interface.aPatchID()].type()
                                    << ", instead of "
                                    << fixedGradientFvPatchScalarField::typeName
                                    << abort(FatalError);
                        }
                    }

                    storage.c().correctBoundaryConditions();
                }
            }

            surfaceScalarField gxf("gxf", g & mesh.Cf());

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

#           include "UpLoop_UEqn.H"

            // --- Pressure corrector PISO loop
            while (control.correct())
            {
#               include "UpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }

#           include "UpLoop_cEqn.H"
#           include "UpLoop_rhoUpdate.H"

            // Make the fluxes absolute
            fvc::makeAbsolute(phi, U);

            interface.updatePoints(); // abs phi

// TODO: Continuity errors
// #           include "trackedSurfaceContinuityErrs.H"
        }
    }
