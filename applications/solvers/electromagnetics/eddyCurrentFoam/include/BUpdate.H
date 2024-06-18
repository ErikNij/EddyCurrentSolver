
// TODO: To save memory, it might be better not to store this field
    // Calculate B
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | BUpdate.H : "
                << "Calculate B."
                << endl;
        }

        globalStorage.BRe()[DEFAULT] = fvc::curl(globalStorage.ARe()[DEFAULT]);
        globalStorage.BIm()[DEFAULT] = fvc::curl(globalStorage.AIm()[DEFAULT]);

        // Update B in conductor region
        globalStorage.BRe().mapInternalField(CONDUCTOR);
        globalStorage.BIm().mapInternalField(CONDUCTOR);

        // If Biot-Savart IS used, A0 in the non-cunducting region
        // is NOT calculated. Thus, we are not allowed to do fvc::curl(A + A0)
        // in the whole default region.
        if (globalSettings.biotSavart)
        {
            globalStorage.BRe()[CONDUCTOR] += fvc::curl(globalStorage.A0Re()[CONDUCTOR]);
            globalStorage.BIm()[CONDUCTOR] += fvc::curl(globalStorage.A0Im()[CONDUCTOR]);
        }

        // Extrapolate boundary field of B
        fvc::extrapolate(globalStorage.BRe()[CONDUCTOR]);
        fvc::extrapolate(globalStorage.BIm()[CONDUCTOR]);
    }
