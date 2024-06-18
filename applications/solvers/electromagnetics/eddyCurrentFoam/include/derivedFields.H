
// TODO: To save memory, it might be better not to store these fields
    // Derived fields in conductor region
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | derivedFields.H : "
                << "Calculate derived fields."
                << endl;
        }

        SM_REGIONSCOPE(CONDUCTOR);

        volScalarField& sigma = globalStorage.sigma()[CONDUCTOR];
        volScalarField& nur = globalStorage.nur()[CONDUCTOR];

        volVectorField& BRe = globalStorage.BRe()[CONDUCTOR];
        volVectorField& BIm = globalStorage.BIm()[CONDUCTOR];

        // Eddy current density
        storage.jRe() =   omega0 * sigma * globalStorage.AIm()[CONDUCTOR];
        storage.jIm() = - omega0 * sigma * globalStorage.ARe()[CONDUCTOR];

        if (control.meshIs3D())
        {
            storage.jRe() -= sigma * globalStorage.VReGrad()[CONDUCTOR];
            storage.jIm() -= sigma * globalStorage.VImGrad()[CONDUCTOR];
        }

        if (globalSettings.biotSavart)
        {
            storage.jRe() += omega0 * sigma * globalStorage.A0Im()[CONDUCTOR];
            storage.jIm() -= omega0 * sigma * globalStorage.A0Re()[CONDUCTOR];
        }

        if (globalStorage.item_F().enabled())
        {
            // Time-averaged Lorentz-force
            globalStorage.F()[CONDUCTOR] = 0.5 * ( (storage.jRe() ^ BRe)
                                                 + (storage.jIm() ^ BIm) );
        }

        if (globalStorage.item_pB().enabled())
        {
            // Time-averaged magnetic pressure
            globalStorage.pB()[CONDUCTOR] = 0.5 * physicalConstant::nu0 * nur
                                          * 0.5 * ( (BRe & BRe)
                                                  + (BIm & BIm) );
        }

        if (globalStorage.item_Q().enabled())
        {
            // Time-averaged joule heating
            globalStorage.Q()[CONDUCTOR] = 0.5 / sigma * ( (storage.jRe() & storage.jRe())
                                                         + (storage.jIm() & storage.jIm()) );
        }
    }

