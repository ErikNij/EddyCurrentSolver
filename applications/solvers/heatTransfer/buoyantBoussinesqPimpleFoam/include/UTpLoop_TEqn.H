
    {
        volScalarField kappaEff("kappaEff", turbulence.nu()/Pr + alphat);

        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(kappaEff, T)
//          ==
//             radiation.ST(rhoCpRef, T)
        );

        if (storage.item_Q().enabled())
        {
            if (Control::debug)
            {
                Info<< Control::typeName << " | TEqn.H : "
                    << "Adding heat source to TEqn."
                    << endl;
            }

            TEqn -= storage.Q()/storage.rhoRef()/storage.CpRef();
        }

        TEqn.relax();

        TEqn.solve();

    //     radiation.correct();
    }
