if (!ddtPhiCorr)
{
    Info << "Correcting interface flux" << endl;

    phi.boundaryField()[interface.aPatchID()] +=
        (
            (
                U.boundaryField()[interface.aPatchID()].snGrad()
               /mesh.boundary()[interface.aPatchID()].deltaCoeffs()
            )
          & mesh.Sf().boundaryField()[interface.aPatchID()]
        );
}
