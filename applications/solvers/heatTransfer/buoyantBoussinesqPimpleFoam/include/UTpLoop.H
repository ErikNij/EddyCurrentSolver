
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | UTpLoop.H : "
                << "Commencing PIMPLE U-T-p loop."
                << endl;
        }

        // --- PIMPLE corrector loop
        while (control.loop())
        {
            uniformDimensionedVectorField& g = storage.g();
            volScalarField& T = storage.T();
            volScalarField& p = storage.p();
            volVectorField& U = storage.U();
            surfaceScalarField& phi = storage.phi();
            uniformDimensionedScalarField& beta = storage.beta();
            uniformDimensionedScalarField& TRef = storage.TRef();
            uniformDimensionedScalarField& Pr = storage.Pr();
            uniformDimensionedScalarField& Prt = storage.Prt();
            incompressible::turbulenceModel& turbulence = storage.turbulence();
            volScalarField& rhok = storage.rhok();
            volScalarField& alphat = storage.alphat();

#           include "UTpLoop_UEqn.H"
#           include "UTpLoop_alphatUpdate.H"
#           include "UTpLoop_TEqn.H"
#           include "UTpLoop_rhokUpdate.H"

            // --- Pressure corrector loop
            while (control.correct())
            {
#               include "UTpLoop_pEqn.H"
            }

            if (control.turbCorr())
            {
                turbulence.correct();
            }
        }
    }
