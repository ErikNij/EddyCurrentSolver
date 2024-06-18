
// TODO: Read only once
    if (storage.item_sigma().enabled())
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | calcSigma.H : "
                << "Calculate elctrical conductivity"
                << endl;
        }

        word phase1Name(storage.transport().phase1Name());
        word phase2Name(storage.transport().phase2Name());

        dimensionedScalar sigma1(storage.transport().subDict(phase1Name).lookup("sigmae"));
        dimensionedScalar sigma2(storage.transport().subDict(phase2Name).lookup("sigmae"));

        volScalarField limitedAlpha1 = min(max(storage.alpha1(), scalar(0)), scalar(1));

        storage.sigma() = limitedAlpha1*sigma1 + (scalar(1) - limitedAlpha1)*sigma2;
        storage.sigma().correctBoundaryConditions();
    }
