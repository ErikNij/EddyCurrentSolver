label nAlphaCorr
(
    readLabel(control.dict().lookup("nAlphaCorr"))
);

label nAlphaSubCycles
(
    readLabel(control.dict().lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaUpLoop_subCycle_alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
#       include "alphaUpLoop_subCycle_alphaEqn.H"
}

interface.correct();

rho == alpha1*transport.rho1() + (scalar(1) - alpha1)*transport.rho2();
