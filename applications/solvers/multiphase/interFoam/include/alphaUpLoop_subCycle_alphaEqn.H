{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);

        rhoPhi = phiAlpha*(transport.rho1() - transport.rho2()) + phi*transport.rho2();
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;

    // Cut compression
// TODO: Compression should depend on time step size
// TODO: How to fix phi after changing aplha?
    dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
    if(pisoDict.lookupOrDefault("cutCompression", false))
    {
        volScalarField magSqrGradAlpha("magSqrGradAlpha", magSqr(fvc::grad(alpha1)));
        magSqrGradAlpha /= max(magSqrGradAlpha)
            + dimensionedScalar(word(), magSqrGradAlpha.dimensions(), VSMALL);

//         volScalarField cutFct("cutFct", 1.0-magSqrGradAlpha);
        volScalarField cutFct("cutFct", 1.0-pow(magSqrGradAlpha, 2));
//         volScalarField cutFct("cutFct", 1.0-pow(magSqrGradAlpha, 4));
        cutFct /= max(cutFct)
            + dimensionedScalar(word(), cutFct.dimensions(), VSMALL);

        volScalarField cutFctAverage("cutFctAverage", fvc::average(fvc::interpolate(cutFct)));

        scalar cAlphaCut = pisoDict.lookupOrDefault("cAlphaCut", 0.0);
        scalar cAlphaCutRelErr =
            pisoDict.lookupOrDefault("cAlphaCutRelErr", 1e-4);
        scalar cAlphaCutPatchEnhance =
            pisoDict.lookupOrDefault("cAlphaCutPatchEnhance", 1.0);

        volScalarField cAlphaCutField(cAlphaCut*cutFctAverage);

        // Shift and symmetrically scale up alpha
        volScalarField newAlpha1
        (
            "newAlpha1",
            (1.0 + cAlphaCutField)*2.0*(alpha1 - 0.5)
        );

        // Patch enhanced cut compression
        const fvBoundaryMesh& boundary = mesh.boundary();
        forAll(boundary, patchi)
        {
            if (isA<wallFvPatch>(boundary[patchi]))
            {
                const fvPatch& patch = boundary[patchi];
                const labelList& faceCells = patch.faceCells();

                forAll(faceCells, facei)
                {
                    label celli = faceCells[facei];

                    scalar cAlphaPatchCut =
                        cAlphaCutPatchEnhance*cAlphaCutField[celli];

                    newAlpha1[celli] =
                        (1.0 + cAlphaPatchCut)*2.0*(alpha1[celli] - 0.5);
                }
            }
        }

        // Cut new alpha
        newAlpha1 = 0.5*min(max(newAlpha1, scalar(-1)), scalar(1)) + 0.5;
        newAlpha1.correctBoundaryConditions();

        // Calculate error in volume fraction of phase1
        dimensionedScalar newAlpha1Vfrac =
            newAlpha1.weightedAverage(mesh.V());
        dimensionedScalar relErrNewAlpha1Vfrac =
            (newAlpha1Vfrac-alpha1Vfrac)/alpha1Vfrac;

        // Correct volume fraction of phase1
        if (mag(relErrNewAlpha1Vfrac).value() > cAlphaCutRelErr)
        {
            newAlpha1 *= (1.0 + cAlphaCutField);

            if (relErrNewAlpha1Vfrac.value() > 0.0)
            {
                newAlpha1 -= cAlphaCutField;
            }

            newAlpha1 = min(max(newAlpha1, scalar(0)), scalar(1));
            newAlpha1.correctBoundaryConditions();

            newAlpha1Vfrac = newAlpha1.weightedAverage(mesh.V());
        }

        alpha1 = newAlpha1;
        alpha1.correctBoundaryConditions();
    }
}
