

    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | alphaUpLoop.H : "
                << "Commencing PIMPLE alpha-U-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            uniformDimensionedVectorField& g = storage.g();
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
            volScalarField& rho = storage.rho();
            surfaceScalarField& rhoPhi = storage.rhoPhi();
            twoPhaseMixture& transport = storage.transport();
            incompressible::turbulenceModel& turbulence = storage.turbulence();
            interfaceProperties& interface = storage.interface();
            volScalarField& alpha1 = storage.alpha1();
            uniformDimensionedScalarField& alpha1Vfrac = storage.alpha1Vfrac();

#           include "alphaUpLoop_subCycle.H"

#           include "alphaUpLoop_UEqn.H"

            // --- Pressure corrector PISO loop
            while (control.correct())
            {
#               include "alphaUpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }

#           include "calcSigma.H"
        }
    }
