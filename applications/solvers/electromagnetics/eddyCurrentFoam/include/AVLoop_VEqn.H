
    // Solve V
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];

        for (int cI = 0; cI < 2; cI++)
        {
            label complexSign = 0;
            volScalarField* VPtr = NULL;
            volVectorField* VGradPtr = NULL;
            volVectorField* APtr = NULL;
            volVectorField* A0Ptr = NULL;
            surfaceScalarField* phiDdtAPtr = NULL;

            if (cI == 0)
            {
                complexSign = -1;
                VPtr = &storage.VRe();
                VGradPtr = &globalStorage.VReGrad()[CONDUCTOR];
                APtr = &globalStorage.AIm()[CONDUCTOR];
                phiDdtAPtr = &storage.phiDdtAIm();

                if (globalSettings.biotSavart)
                {
                    A0Ptr = &globalStorage.A0Im()[CONDUCTOR];
                }
            }
            else
            {
                complexSign = 1;
                VPtr = &storage.VIm();
                VGradPtr = &globalStorage.VImGrad()[CONDUCTOR];
                APtr = &globalStorage.ARe()[CONDUCTOR];
                phiDdtAPtr = &storage.phiDdtARe();

                if (globalSettings.biotSavart)
                {
                    A0Ptr = &globalStorage.A0Re()[CONDUCTOR];
                }
            }

            volScalarField& V = *VPtr;
            volVectorField& VGrad = *VGradPtr;
            volVectorField& A = *APtr;
            surfaceScalarField& phiDdtA = *phiDdtAPtr;

// TODO: Use dedicated boundary condition?
// TODO: Remove interfaceLabel? Use dedicated boundary condition!
            // Update gradient boundary condition for V
            {
                label patchI = control.interfaceLabel();

                if (isA<fixedGradientFvPatchScalarField>(V.boundaryField()[patchI]))
                {
                    fixedGradientFvPatchScalarField& VP =
                        refCast<fixedGradientFvPatchScalarField>
                        (
                            V.boundaryField()[patchI]
                        );

                    // Set correct gradient
                    VP.gradient() = -phiDdtA.boundaryField()[patchI]/VP.patch().magSf();
                }

                V.correctBoundaryConditions();
            }

            surfaceScalarField sigmaf("sigmaf", fvc::interpolate(sigma));

            // Read Laplacian scheme for discretising laplacian(sigmaf,V)
            word laplacianSchemeV
            (
                mesh.schemesDict().laplacianScheme
                (
                    "laplacian(" + sigmaf.name() + "," + V.name() + ")"
                )
            );

            // Prepare V-equation
            fvScalarMatrix VEqn(fvm::laplacian(sigmaf, V));

            if (laplacianSchemeV != "jumpGauss")
            {
                // NOTE: We do NOT use sigmaf here!
                volVectorField sigmaGrad = fvc::grad(sigma);

                // Explicit source term for V due to variable sigma
                VEqn += complexSign * omega0 * (sigmaGrad & A);

                if (globalSettings.biotSavart)
                {
                    volVectorField& A0 = *A0Ptr;

                    VEqn += complexSign * omega0 * (sigmaGrad & A0);
                }
            }

            // Read solver for V
            word solverV
            (
                mesh.solutionDict().solver(V.name()).lookup("solver")
            );

            if(solverV != "deflatedPCG")
            {
// TODO: Read reference values
//                 VEqn.setReference(VRefCell, VRefValue);
                VEqn.setReference(0, 0.0);
            }

            // Solve V-equation
            VEqn.solve(control.subDict(Region::CONDUCTOR, V.name()));

            if (laplacianSchemeV != "jumpGauss")
            {
                VGrad = fvc::reconstruct(VEqn.flux()/sigmaf);
            }
            else
            {
                VGrad =
                    fvc::reconstruct
                    (
                        VEqn.flux()/ownerInterpolate(sigma),
                        VEqn.jumpFlux()/neighbourInterpolate(sigma)
                    );
            }
        }
    }
