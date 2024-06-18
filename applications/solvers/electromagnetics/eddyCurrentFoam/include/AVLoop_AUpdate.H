
    // Update flux of gradient of normal component of A
    // in default region
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(DEFAULT);

        for (int cI = 0; cI < 2; cI++)
        {
            volVectorField* APtr = NULL;
            surfaceVectorField* phiGradAnPtr = NULL;

            if (cI == 0)
            {
                APtr = &globalStorage.ARe()[DEFAULT];
                phiGradAnPtr = &storage.phiGradAnRe();
            }
            else
            {
                APtr = &globalStorage.AIm()[DEFAULT];
                phiGradAnPtr = &storage.phiGradAnIm();
            }

            volVectorField& A = *APtr;
            surfaceVectorField& phiGradAn = *phiGradAnPtr;

            surfaceTensorField gradAf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A)),
                    "interpolate(grad(" + A.name() + "))"
                );

            phiGradAn = -gradAf & mesh.Sf();

            phiGradAn += tr(gradAf) * mesh.Sf();
        }
    }

    // Update A in conductor region
    // from data in default region
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | AUpdate.H : "
                << "Update A in conductor region."
                << endl;
        }

// TEST: Use relaxation?
//         globalStorage.ARe()[CONDUCTOR].storePrevIter();
//         globalStorage.ARe()[CONDUCTOR].storePrevIter();

        // Update A
        globalStorage.ARe().map(CONDUCTOR);
        globalStorage.AIm().map(CONDUCTOR);

// TEST: Use relaxation?
//         const scalar AReRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.ARe()[DEFAULT].name());
//         const scalar AImRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.AIm()[DEFAULT].name());
//         globalStorage.ARe()[CONDUCTOR].relax(AReRf);
//         globalStorage.ARe()[CONDUCTOR].relax(AImRf);
    }


    // Update flux of rate of change of A
    // in conductor region
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

// TEST/TODO: Interpolation with nuf?
// TEST/TODO: Better map phiGradAn from default region instead of re-calculation?
        for (int cI = 0; cI < 2; cI++)
        {
            word  complexPart("");
            label complexSign = 0;
            volVectorField* APtr = NULL;
            surfaceScalarField* phiDdtAPtr = NULL;

            if (cI == 0)
            {
                complexPart = "Re";
                complexSign = 1;
                APtr = &globalStorage.ARe()[CONDUCTOR];
                phiDdtAPtr = &storage.phiDdtARe();
            }
            else
            {
                complexPart = "Im";
                complexSign = -1;
                APtr = &globalStorage.AIm()[CONDUCTOR];
                phiDdtAPtr = &storage.phiDdtAIm();
            }

            volVectorField& A = *APtr;
            surfaceScalarField& phiDdtA = *phiDdtAPtr;

            surfaceTensorField gradAf =
                fvc::interpolate
                (
                    fvc::grad(linearInterpolate(A)),
                    "interpolate(grad(" + A.name() + "))"
                );

            surfaceVectorField phiGradAn("phiGradAn" + complexPart, -gradAf & mesh.Sf());

            phiGradAn += tr(gradAf) * mesh.Sf();

            phiDdtA = complexSign * omega0 * (mesh.Sf() & fvc::interpolate(A));
        }

// TEST/TODO: Interpolation with nuf?
//         storage.phiDdtARe() =  omega0*(mesh.Sf() & fvc::interpolate(globalStorage.ARe()[CONDUCTOR]));
//         storage.phiDdtAIm() = -omega0*(mesh.Sf() & fvc::interpolate(globalStorage.AIm()[CONDUCTOR]));

        if (globalSettings.biotSavart)
        {
            storage.phiDdtARe() += omega0*(mesh.Sf() & fvc::interpolate(globalStorage.A0Re()[CONDUCTOR]));
            storage.phiDdtAIm() -= omega0*(mesh.Sf() & fvc::interpolate(globalStorage.A0Im()[CONDUCTOR]));
        }
    }
