
    // Solve A0 using Biot-Savart law
    if (globalSettings.biotSavart)
    {
        // Calculate impressed vector potential A0
        {
            SM_REGIONSCOPE(CONDUCTOR);

            if (Control::debug)
            {
                Info<< Control::typeName << " | A0BiotSavart.H : "
                    << "Calculate impressed vector potential"
                    << " via Biot-Savart's law."
                    << endl;
            }

// TODO: Put this into manager! and write only once
            // Construct Biot-Savart class
            edgeBiotSavart eBiotSavart(globalMesh[CONDUCTOR]);
//             eBiotSavart.writeVTK();

            // Calculate A0
            eBiotSavart.A(globalStorage.A0Re()[CONDUCTOR], edgeBiotSavart::REAL);
            eBiotSavart.A(globalStorage.A0Im()[CONDUCTOR], edgeBiotSavart::IMAG);
        }

        // Update A0 in default region
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | A0BiotSavart.H : "
                    << "Update A0 in default region."
                    << endl;
            }

            // Update A0
            globalStorage.A0Re().rmap(CONDUCTOR);
            globalStorage.A0Im().rmap(CONDUCTOR);
        }

        // Calculate effective source current for magnetization
        if (control.meshIs3D())
        {
            SM_REGIONSCOPE(CONDUCTOR);

            volScalarField chiNur
            (
                IOobject
                (
                    "chiNur",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(word(), dimless, SMALL)
            );

            chiNur.internalField() =
                1.0 - globalStorage.nur()[CONDUCTOR].internalField() + SMALL;

            // Correct boundary conditions
            // NOTE: Important for parallel runs
            chiNur.correctBoundaryConditions();

            surfaceScalarField chiNurf("chiNurf", fvc::interpolate(chiNur));

            // Calculate flux of gradient of normal component of A0
            // and effective source current
            for (int cI = 0; cI < 2; cI++)
            {
                word  complexPart("");
                volVectorField* A0Ptr = NULL;
                regionVolVectorField* j0Ptr = NULL;

                if (cI == 0)
                {
                    complexPart = "Re";
                    A0Ptr = &globalStorage.A0Re()[CONDUCTOR];
                    j0Ptr = &globalStorage.j0Re();
                }
                else
                {
                    complexPart = "Im";
                    A0Ptr = &globalStorage.A0Im()[CONDUCTOR];
                    j0Ptr = &globalStorage.j0Im();
                }

                volVectorField& AM = *A0Ptr;
                regionVolVectorField& jM = *j0Ptr;

                word A0name = AM.name();
                AM.rename("AM" + complexPart);

                surfaceTensorField gradAMf =
                    fvc::interpolate
                    (
                        fvc::grad(linearInterpolate(AM)),
                        "interpolate(grad(" + AM.name() + "))"
                    );

                surfaceVectorField phiGradAMn
                (
                    "phiGradAMn" + complexPart,
                    -gradAMf & mesh.Sf()
                );

                phiGradAMn += tr(gradAMf) * mesh.Sf();

                // Make sure the source current is zero in default region, first
                jM[DEFAULT] = dimensionedVector(word(), dimCurrent/dimArea, vector::zero);

// TEST:
// fvc::curl(chiNur*(2.0*(*skew(fvc::grad(AM, word("curl(" + A0name + ")"))))))().write();
// fvc::laplacian(chiNurf, AM)().write();

                // Effective source current for magnetization
                jM[CONDUCTOR] = - physicalConstant::nu0 * fvc::laplacian(chiNurf, AM);

// TEST:
// jM[CONDUCTOR].write();

                AM.rename(A0name);
            }

            // Update j0 in default region
            {
                if (Control::debug)
                {
                    Info<< Control::typeName << " | A0BiotSavart.H : "
                        << "Update j0 in default region."
                        << endl;
                }

                // Update j0
                globalStorage.j0Re().rmap(CONDUCTOR);
                globalStorage.j0Im().rmap(CONDUCTOR);
            }
        }
    }
