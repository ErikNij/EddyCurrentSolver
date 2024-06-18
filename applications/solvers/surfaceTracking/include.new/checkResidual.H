{
    // const scalarField& pI = p.internalField();
    // const scalarField& pOI = p.prevIter().internalField();

    dimensionedScalar residual = 
        (gMax((p - p.prevIter())()) - gMin((p - p.prevIter())()))
       /(gMax(p) - gMin(p) + SMALL);

    Info << "Pressure residual: " << residual.value() 
        << ", nIter: " << corr << endl;

    if (residual.value() < 1e-12)
    {
         break;
    }
}
