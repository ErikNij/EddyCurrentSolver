
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UpLoop.H : "
                << "Commencing PIMPLE U-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
            incompressible::turbulenceModel& turbulence = storage.turbulence();

            #include "UpLoop_UEqn.H"

            // --- Pressure corrector loop
            while (control.correct())
            {
                #include "UpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }
        }
    }
