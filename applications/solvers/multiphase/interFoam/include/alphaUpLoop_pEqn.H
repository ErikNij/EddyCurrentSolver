
    {
        volScalarField rAU("rAU", 1.0/UEqn.A());
        surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

        U = rAU*UEqn.H();

        surfaceScalarField phiU
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

        phi = phiU +
            (
                fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
              - (g & mesh.Cf())*fvc::snGrad(rho)
            )*rAUf*mesh.magSf();

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
                            << "is being corrected by means of a "
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

// TODO: Test
        // Stop phase2 completely
        dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
        if(pisoDict.lookupOrDefault("phase1Only", false))
        {
            scalar phase1Cut = pisoDict.lookupOrDefault("phase1Cut", 0.01);
            scalar phase1CutOff = pisoDict.lookupOrDefault("phase1CutOff", 0.1*phase1Cut);

            const unallocLabelList& owner = mesh.owner();
            const unallocLabelList& neighbour = mesh.neighbour();

            // Calculate internal correction
            forAll (owner, facei)
            {
                label own = owner[facei];
                label nei = neighbour[facei];

                if(alpha1[own] < phase1Cut && alpha1[nei] < phase1Cut)
                {
                    U[own] = vector::zero;
                    phi[facei] = scalar(0);
                    rhoPhi[facei] = scalar(0);
                }
            }
            forAll(mesh.C(), celli)
            {
                if(alpha1[celli] < phase1CutOff)
                {
                    U[celli] = vector::zero;
                }
            }
            U.correctBoundaryConditions();
        }

// TODO: Continuity errors
// #       include "continuityErrs.H"
    }

