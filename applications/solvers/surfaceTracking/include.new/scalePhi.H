// Info << interface.twoFluids() << ", " << pNeedRef << endl;

if(interface.twoFluids() && pNeedRef)
{
    phi.boundaryField()[interface.bPatchID()] =
        (
            U.boundaryField()[interface.bPatchID()] 
          & mesh.Sf().boundaryField()[interface.bPatchID()]
        );
        
    scalarField weights = 
        mag(phi.boundaryField()[interface.bPatchID()]);

    if(mag(gSum(weights)) > VSMALL)
    {
        weights /= gSum(weights);
    }

    phi.boundaryField()[interface.bPatchID()] -=
        weights*gSum(phi.boundaryField()[interface.bPatchID()]);

    Info << "Free surface shadow flow rate: "
        << gSum(phi.boundaryField()[interface.bPatchID()]) << endl;

    phi.boundaryField()[interface.bPatchID()] +=
        p.boundaryField()[interface.bPatchID()].snGrad()
       *mesh.magSf().boundaryField()[interface.bPatchID()]
       /AU.boundaryField()[interface.bPatchID()];
}


// Scale phi for fluid A patches if necessary
if (false)
{
    // Check if fluid A needs reference
    bool needRef = true;
    forAll (p.boundaryField(), patchI)
    {
        // Check only fluid A
        if (patchI != interface.aPatchID())
        {
            if (fluidIndicator.boundaryField()[patchI].size())
            {
                if (min(fluidIndicator.boundaryField()[patchI]) > SMALL)
                {
                    if (p.boundaryField()[patchI].fixesValue())
                    {
                        needRef = false;
                    }
                }
            }
        }
    }
    needRef = returnReduce<bool>(needRef, andOp<bool>());

    if (needRef)
    {
        Info << "Scaling phi for fluid A" << endl;

        scalar massIn = 0.0;
        scalar fixedMassOut = 0.0;
        scalar adjustableMassOut = 0.0;

        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            const fvsPatchScalarField& phip = phi.boundaryField()[patchi];

            if 
            (
                !Up.coupled() 
             && (patchi != interface.aPatchID())
            )
            {
                bool fluidAPatch = false;

                if (fluidIndicator.boundaryField()[patchi].size())
                {
                    if 
                    (
                        min(fluidIndicator.boundaryField()[patchi]) 
                      > SMALL
                    )
                    {
                        fluidAPatch = true;
                    }
                }

                fluidAPatch = returnReduce<bool>(fluidAPatch, orOp<bool>());
                
                if (fluidAPatch)
                {
//                     Info << "Fluid A patch: " 
//                         << mesh.boundary()[patchi].name() << endl;
                    if
                    (
                        Up.fixesValue()
                    )
                    {
                        forAll (phip, i)
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                fixedMassOut += phip[i];
                            }
                        }
                    }
                    else
                    {
                        forAll(phip, i)
                        {
                            if (phip[i] < 0.0)
                            {
                                massIn -= phip[i];
                            }
                            else
                            {
                                adjustableMassOut += phip[i];
                            }
                        }
                    }
                }
            }
        }

        reduce(massIn, sumOp<scalar>());
        reduce(fixedMassOut, sumOp<scalar>());
        reduce(adjustableMassOut, sumOp<scalar>());

        scalar massCorr = 1.0;

        if (mag(adjustableMassOut) > SMALL)
        {
            massCorr = (massIn - fixedMassOut)/adjustableMassOut;
        }

        Info<< "bool Foam::adjustPhi(...) massIn: " << massIn
            << " fixedMassOut: " << fixedMassOut
            << " adjustableMassOut: " << adjustableMassOut
            << " mass corr: " << massCorr
            << endl;
        
        forAll (phi.boundaryField(), patchi)
        {
            const fvPatchVectorField& Up = U.boundaryField()[patchi];
            fvsPatchScalarField& phip = phi.boundaryField()[patchi];

            if (!Up.coupled() && (patchi != interface.aPatchID()))
            {
                bool fluidAPatch = false;

                if (fluidIndicator.boundaryField()[patchi].size())
                {
                    if 
                    (
                        min(fluidIndicator.boundaryField()[patchi]) 
                      > SMALL
                    )
                    {
                        fluidAPatch = true;
                    }
                }

                fluidAPatch = returnReduce<bool>(fluidAPatch, orOp<bool>());
                
                if (fluidAPatch)
                {
                    if
                    (
                        !Up.fixesValue()
                    )
                    {
                        forAll (phip, i)
                        {
                            if (phip[i] > 0.0)
                            {
                                phip[i] *= massCorr;
                            }
                        }
                    }
                }
            }
        }
    }
}
