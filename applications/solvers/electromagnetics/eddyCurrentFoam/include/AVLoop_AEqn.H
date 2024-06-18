
    // Solve A component-wise in default region
    if (!(globalSettings.biotSavart && globalSettings.lowFrequency))
    {
        SM_REGIONSCOPE(DEFAULT);

        volScalarField& sigma = globalStorage.sigma()[DEFAULT];
        volScalarField& nur = globalStorage.nur()[DEFAULT];

        volVectorField& ARe = globalStorage.ARe()[DEFAULT];
        volVectorField& AIm = globalStorage.AIm()[DEFAULT];

        volVectorField& j0Re = globalStorage.j0Re()[DEFAULT];
        volVectorField& j0Im = globalStorage.j0Im()[DEFAULT];

        // Interpolate relative magnetic reluctivity
        surfaceScalarField nurf("nurf", fvc::interpolate(nur));

        // Assemble base-matrices for A
        fvVectorMatrix AReEqn
        (
          - fvm::laplacian(nurf,ARe)
        );

        fvVectorMatrix AImEqn
        (
          - fvm::laplacian(nurf,AIm)
        );

        // Add source current
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add source current density to AEqn"
                    << endl;
            }

            AReEqn -= physicalConstant::mu0 * j0Re;
            AImEqn -= physicalConstant::mu0 * j0Im;
        }

        // Add source vector potential if Biot-Savart is used
        if (globalSettings.biotSavart)
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add impressed vector potenial to AEqn"
                    << endl;
            }

            AReEqn -= physicalConstant::mu0 * sigma * omega0 * globalStorage.A0Im()[DEFAULT];
            AImEqn += physicalConstant::mu0 * sigma * omega0 * globalStorage.A0Re()[DEFAULT];
        }

        // Add explicit gradient terms of the
        // electric scalar potential V
        if (control.meshIs3D() && !globalSettings.lowFrequency)
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | AEqn.H : "
                    << "Add gradient of electric scalar potential to AEqn"
                    << endl;
            }

            AReEqn += physicalConstant::mu0 * sigma * globalStorage.VReGrad()[DEFAULT];
            AImEqn += physicalConstant::mu0 * sigma * globalStorage.VImGrad()[DEFAULT];
        }

// TODO: Which discretization scheme is used? If not jumpLaplacian, we need to add gradient terms!

        // Create list with boundary types
        wordList ARePatchTypes(mesh.boundary().size());
        wordList AImPatchTypes(mesh.boundary().size());

        forAll(ARe.boundaryField(), patchI)
        {
            ARePatchTypes[patchI] = ARe.boundaryField()[patchI].type();
        }
        forAll(AIm.boundaryField(), patchI)
        {
            AImPatchTypes[patchI] = AIm.boundaryField()[patchI].type();
        }

        // Loop over components
        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            // Only solve for valid directions
            if(control.solutionDir()[cmpt] < 1) continue;

            // Create component fields
            volScalarField AReCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]) + "Re",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                ARe.component(cmpt),
                ARePatchTypes
            );
            volScalarField AImCmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]) + "Im",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                AIm.component(cmpt),
                AImPatchTypes
            );

            // Correct/Evaluate boundary values
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Create temporary component equations
            tmp<fvScalarMatrix> tAReCmptEqn
            (
                new fvScalarMatrix
                (
                    AReCmpt,
                    AReEqn.dimensions()
                )
            );
            tmp<fvScalarMatrix> tAImCmptEqn
            (
                new fvScalarMatrix
                (
                    AImCmpt,
                    AImEqn.dimensions()
                )
            );

            // Get refs to component equations
            fvScalarMatrix& AReCmptEqn = tAReCmptEqn();
            fvScalarMatrix& AImCmptEqn = tAImCmptEqn();

            // Copy symmetric matrix coeffs and source
            AReCmptEqn.diag() = AReEqn.diag();
            AReCmptEqn.upper() = AReEqn.upper();
            AReCmptEqn.internalCoeffs() = AReEqn.internalCoeffs().component(cmpt);
            AReCmptEqn.boundaryCoeffs() = AReEqn.boundaryCoeffs().component(cmpt);
            AReCmptEqn.source() = AReEqn.source().component(cmpt);
            AImCmptEqn.diag() = AImEqn.diag();
            AImCmptEqn.upper() = AImEqn.upper();
            AImCmptEqn.internalCoeffs() = AImEqn.internalCoeffs().component(cmpt);
            AImCmptEqn.boundaryCoeffs() = AImEqn.boundaryCoeffs().component(cmpt);
            AImCmptEqn.source() = AImEqn.source().component(cmpt);

            // Create complex vector2 field
            volVector2Field ACmpt
            (
                IOobject
                (
                    "A" + word(vector::componentNames[cmpt]),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector2(word(), dimless, vector2::zero)
            );

            // Init complex vector2 field
            ACmpt.replace(0, AReCmpt);
            ACmpt.replace(1, AImCmpt);

            // Prepare block system for current component
            fvBlockMatrix<vector2> AcmptEqn(ACmpt);

            // Lookup solver type for AcmptEqn
            word solverType =
                mesh.solutionDict().subDict("solvers").
                subDict(AcmptEqn.psi().name()).lookup("solver");

            // Force diag/upper/lower square for AMG (matrix for A is asymmetric)
            if (solverType == "AMG")
            {
                Field<tensor2>& AEqnD = AcmptEqn.diag().asSquare(); AEqnD = tensor2::zero;
                Field<tensor2>& AEqnU = AcmptEqn.upper().asSquare(); AEqnU = tensor2::zero;
                Field<tensor2>& AEqnL = AcmptEqn.lower().asSquare(); AEqnL = tensor2::zero;
            }

            // Insert equations into block matrix
            AcmptEqn.insertEquation(0, AReCmptEqn);
            AcmptEqn.insertEquation(1, AImCmptEqn);

            // Clear memory
            tAReCmptEqn.clear();
            tAImCmptEqn.clear();

            // Add coupling terms for real/imagninary ARe-AIm
            // dependence of magnetic vector potential
            if (!globalSettings.lowFrequency)
            {
                AcmptEqn.insertEquationCoupling(0, 1, - physicalConstant::mu0 * sigma * omega0);
                AcmptEqn.insertEquationCoupling(1, 0,   physicalConstant::mu0 * sigma * omega0);
            }

            // Solve A but prevent segfault for homogenous solution
            scalar preResidual = mag(gSum(cmptMag(AcmptEqn.residual(ACmpt, AcmptEqn.source()))));
            if (preResidual > 0)
            {
                AcmptEqn.solve(control.subDict(Region::DEFAULT, ACmpt.name()));
            }

            // Retrieve solution
            AcmptEqn.retrieveSolution(0, AReCmpt.internalField());
            AcmptEqn.retrieveSolution(1, AImCmpt.internalField());

            // Correct component boundary conditions
            AReCmpt.correctBoundaryConditions();
            AImCmpt.correctBoundaryConditions();

            // Replace components of A
            ARe.replace(cmpt, AReCmpt);
            AIm.replace(cmpt, AImCmpt);
        }

        // Correct boundary conditions
        ARe.correctBoundaryConditions();
        AIm.correctBoundaryConditions();
    }
