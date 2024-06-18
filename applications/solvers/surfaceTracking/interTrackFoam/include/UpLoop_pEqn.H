
    {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        U = rAU*UEqn.H();

        surfaceScalarField phiU
        (
            (fvc::interpolate(U) & mesh.Sf())
// TODO: Port ddtCorr from OpenFOAM 3.0.x or use new approach from Hrv
//           + fvc::ddtPhiCorr(rAU, U, phi)
        );

        phi = phiU;

        phi -= rAUf*gxf*fvc::snGrad(rho)*mesh.magSf();

        if (storage.item_F().enabled())
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | pEqn.H : "
                    << "Adding volume force flux from F to phi."
                    << endl;
            }

            phi += rAUf*(fvc::interpolate(storage.F()) & mesh.Sf());
        }

// TODO: Two fluids
// #       include "scalePhi.H"
//     if (p.needReference())
//     {
//         fvc::makeRelative(phiHbyA, U);
//         adjustPhi(phiHbyA, U, p);
//         fvc::makeAbsolute(phiHbyA, U);
//     }

        // Manually correct the boundary conditions for the pressure at walls
// TODO: Do not use fixedGradient
        forAll(p.boundaryField(), patchI)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchI]))
            {
                if
                (
                    p.boundaryField()[patchI].type()
                    == fixedFluxPressureFvPatchScalarField::typeName
                )
                {
                    if (Control::debug)
                    {
                        Info<< Control::typeName << " | pEqn.H : "
                            << "Pressure at patch "
                            << mesh.boundary()[patchI].name()
                            << " is being corrected by means of a "
                            << "'fixedFluxPressure' boundary condition. "
                            << "This is OK."
                            << endl;
                    }
                }
                else if
                (
                    p.boundaryField()[patchI].type()
                    == fixedGradientFvPatchScalarField::typeName
                )
                {
                    fixedGradientFvPatchScalarField& pp =
                        refCast<fixedGradientFvPatchScalarField>
                        (
                            p.boundaryField()[patchI]
                        );

                    // Patch data
                    const vectorField& Sfp =
                        mesh.Sf().boundaryField()[patchI];
                    const scalarField& magSfp =
                        mesh.magSf().boundaryField()[patchI];


                    if (!settings.snGradpFromFlux)
                    {
                        if (Control::debug)
                        {
                            Info<< Control::typeName << " | pEqn.H : "
                                << "Correcting pressure at patch "
                                << mesh.boundary()[patchI].name() << "."
                                << endl;
                        }

                            // Direct way by means of F
                        if (storage.item_F().enabled())
                        {
                            const vectorField& Fp =
                                storage.F().boundaryField()[patchI];

                            pp.gradient() = (Fp & Sfp)/magSfp;
                        }
                        else
                        {
                            pp.gradient() = 0.0;
                        }
                    }
                    else
                    {
                        if (Control::debug)
                        {
                            Info<< Control::typeName << " | pEqn.H : "
                                << "Correcting pressure at patch "
                                << mesh.boundary()[patchI].name()
                                << " from flux-difference."
                                << endl;
                        }

                        // Indirect way by means of fixed flux-pressure
                        const scalarField& phip =
                            phi.boundaryField()[patchI];
                        const vectorField& Up =
                            U.boundaryField()[patchI];
                        const scalarField& rAUfp =
                            rAUf.boundaryField()[patchI];

                        pp.gradient() = ((phip - (Up & Sfp))/rAUfp)/magSfp;
                    }

                    pp.updateCoeffs();
                }
                else
                {
                    if (storage.item_F().enabled())
                    {
                        WarningIn(Control::typeName + " | pEqn.H")
                            << "The pressure boundary condition for "
                            << "wall-patch named '"
                            << mesh.boundary()[patchI].name() << "' is '"
                            << p.boundaryField()[patchI].type() << "',"
                            << endl << "    "
                            << "It should be of type 'fixedGradient' "
                            << "or 'fixedFluxPressure' "
                            << "to correctly represent the gradient dp/dn."
                            << endl;
                    }
                }
            }
        }

        p.correctBoundaryConditions();

        // Non-orthogonal pressure corrector loop
        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, p) == fvc::div(phi)
            );

// TODO: Two fluids: Pressure equation reference
// #           include "setReference.H"
//             pEqn.setReference(pRefCell, pRefValue);
            if (!interface.twoFluids())
            {
                // Read solver for p
                word pSolver
                (
                    mesh.solutionDict().solver(p.name()).lookup("solver")
                );

                if(pSolver != "deflatedPCG")
                {
                    pEqn.setReference(control.pRefCell(), control.pRefValue());
                }
            }

            pEqn.solve(mesh.solutionDict().solver(p.select(control.finalInnerIter())));

            // Manually store residuls after first
            // non-orthogonal iteration
            if (control.firstInnerIter())
            {
                control.storeResiduals(p.name());
            }

            // Stop non-orthogonal correction after the next (final)
            // iteration, if the last number of solution iterations
            // was "0" for the named field
            control.skipZeroNonOrtho(p.name());

            if (control.finalNonOrthogonalIter())
            {
                // Store old phi
                surfaceScalarField phiRel("phiRel", phi);

                // Correct phi
                phi -= pEqn.flux();

                // Explicitly relax p and corresponding relaxed phi
                if (!control.finalIter())
                {
                    p.relax();
                }
                phiRel -= pEqn.flux();

                // Momentum corrector with relaxed phi
                U += rAU*fvc::reconstruct((phiRel - phiU)/rAUf);
                U.correctBoundaryConditions();
            }
        }

// TODO: Continuity errors
// #       include "continuityErrs.H"

// TODO: Port ddtCorr from OpenFOAM 3.0.x or use new approach from Hrv
//         {
//             Uf = fvc::interpolate(U);
//             surfaceVectorField n(mesh.Sf()/mesh.magSf());
//             Uf += n*(phi/mesh.magSf() - (n & Uf));
//         }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);
    }

