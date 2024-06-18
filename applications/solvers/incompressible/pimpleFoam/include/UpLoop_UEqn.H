    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence.divDevReff()
    );

    UEqn.relax();

    if (control.momentumPredictor())
    {
        if (storage.item_F().enabled())
        {
            solve
            (
                UEqn == fvc::reconstruct
                (
                  - (fvc::snGrad(p)*mesh.magSf())
                  + (fvc::interpolate(storage.F())/storage.rho() & mesh.Sf())
                )
            );
        }
        else
        {
            solve
            (
                UEqn == fvc::reconstruct
                (
                  - (fvc::snGrad(p)*mesh.magSf())
                )
            );
        }
    }
