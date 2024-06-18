
    scalar CoNum = max
    (
         mesh.surfaceInterpolation::deltaCoeffs()
        *mag(phiNet/fvc::interpolate(rho))/mesh.magSf()
    ).value()*runTime.deltaT().value();

    Info<< "\nMax Courant Number = " << CoNum << endl << endl;

