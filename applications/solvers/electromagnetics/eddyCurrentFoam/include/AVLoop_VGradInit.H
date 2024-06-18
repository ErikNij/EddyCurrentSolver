
    // Update gradient of V in default region
    // from data in conductor region
    if (control.meshIs3D())
    {
        SM_REGIONSCOPE(CONDUCTOR);

        if (Control::debug)
        {
            Info<< Control::typeName << " | VGradInit.H : "
                << "Init gradient of V."
                << endl;
        }

        // Init gradient of V
        globalStorage.VReGrad()[CONDUCTOR] = fvc::grad(storage.VRe());
        globalStorage.VImGrad()[CONDUCTOR] = fvc::grad(storage.VIm());
    }
