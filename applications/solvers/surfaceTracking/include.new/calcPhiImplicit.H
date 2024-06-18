{
    phi = (fvc::interpolate(HU)/fvc::interpolate(AU)) & mesh.Sf();
    
    forAll(phi.boundaryField(), patchI)
    {
        if (!phi.boundaryField()[patchI].coupled())
        {
            phi.boundaryField()[patchI] =
            (
                U.boundaryField()[patchI] 
              & mesh.Sf().boundaryField()[patchI]
            );
        }
        else if (isA<twoPhaseInterfaceVelocityFvPatchVectorField>(U))
        {
            phi.boundaryField()[patchI] =
            (    
                HU.boundaryField()[patchI].patchInternalField()
               /AU.boundaryField()[patchI].patchInternalField()
            ) & mesh.Sf().boundaryField()[patchI];
        }
    }

    if (ddtPhiCorr)
    {
        phi += fvc::ddtPhiCorr(1.0/AU, rho, U, phi);
    }

    // Modify rAU at the interface
    forAll(rAU.boundaryField(), patchI)
    {
        // Calculate rAU at the interface using harmonic
        // mean interpolation
        if
        (
            isA<twoPhaseInterfaceVelocityFvPatchVectorField>
            (
                U.boundaryField()[patchI]
            )
        )
        {
            const scalarField& w = 
                mesh.boundary()[patchI].weights();

            rAU.boundaryField()[patchI] = 
                1.0
               /(
                   w*AU.boundaryField()[patchI].patchNeighbourField()
                 + (1.0 - w)*AU.boundaryField()[patchI].patchInternalField()
                );
        }
    }   
}
