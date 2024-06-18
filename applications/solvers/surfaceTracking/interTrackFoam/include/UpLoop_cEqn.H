
    if (storage.item_c().enabled())
    {
        const uniformDimensionedScalarField& Dc = storage.Dc();
        volScalarField& c = storage.c();

        fvScalarMatrix cEqn
        (
            fvm::ddt(c)
          + fvm::div(phi, c)
          - fvm::laplacian(Dc, c)
        );

        if (settings.relToUinf)
        {
            cEqn += (U-storage.Uinf()) & storage.gradcInf();
        }

        cEqn.relax();

        cEqn.solve();
    }
