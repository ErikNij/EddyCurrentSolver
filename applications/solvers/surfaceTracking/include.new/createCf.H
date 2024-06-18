
// Reading/calculating Cf field for ddtPhiCorr
IOobject CfHeader
(
    "Cf",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

surfaceVectorField Cf
(
    CfHeader,
    mesh,
    dimensionedVector("0", dimLength, vector::zero)
);

if(!CfHeader.headerOk())
{
#   include "updateCf.H"
}
