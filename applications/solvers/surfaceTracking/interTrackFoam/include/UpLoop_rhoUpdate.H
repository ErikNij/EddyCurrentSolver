

// TODO: Two fluids: Use fluidIndicator!
    rho = transport.rho();

    if (storage.item_c().enabled())
    {
        volScalarField& c = storage.c();

        const dictionary& concentrationDict =
            transport.subDict("concentration");

        dimensionedScalar ddcRho(concentrationDict.lookup("ddcRho"));
        dimensionedScalar cRefRho(concentrationDict.lookup("cRefRho"));

// TODO: Two fluids: Use fluidIndicator!
        rho += ddcRho*(c - cRefRho);
    }
