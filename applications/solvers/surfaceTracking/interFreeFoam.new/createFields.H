
Info<< "\nReading field p\n" << endl << flush;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl << flush;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#   include "createPhi.H"
phi.oldTime();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimMass/dimVolume, 0)
);

freeSurface interface(mesh, rho, U, p, phi);

volScalarField* fluidIndicatorPtr = NULL;

if(args.options().found("parallel"))
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}
else
{
    fluidIndicatorPtr = new volScalarField
    (
        IOobject
        (
            "fluidIndicator",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        interface.fluidIndicator()
    );
    fluidIndicatorPtr->write();
}

volScalarField& fluidIndicator = *fluidIndicatorPtr;
fluidIndicator.correctBoundaryConditions();



rho = fluidIndicator
   *(
        interface.rhoFluidA() 
      - interface.rhoFluidB()
    ) 
  + interface.rhoFluidB();
rho.correctBoundaryConditions();

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidIndicator
   *(
        interface.muFluidA() 
      - interface.muFluidB()
    ) 
  + interface.muFluidB()
);
mu.correctBoundaryConditions();

// volScalarField nu("nu", mu/rho);
// nu.correctBoundaryConditions();

dimensionedScalar nu = interface.muFluidA()/interface.rhoFluidA();
nu = dimensionedScalar("nu", nu.dimensions(), nu.value());

#   include "createSf.H"

// Read interpolators if present
if (interface.twoFluids())
{
    interface.interpolatorAB();
}


#   include "setRefCell.H"


areaScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    interface.aMesh().faceCurvatures()
);
K.write();
// const scalarField& K = interface.aMesh().faceCurvatures().internalField();

// const scalarField X = 
//     interface.aMesh().areaCentres().internalField().component(0);

Info << "Free surface curvature: min = " << gMin(K)
    << ", max = " << gMax(K)
    << ", average = " << gAverage(K) << endl << flush;

// forAll(K, faceI)
// {
//     Pout << X[faceI] << " " << K[faceI] << endl;
// }

// return(0);


// volScalarField dist
// (
//     IOobject
//     (
//         "dist",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("0", dimless, 0)
// );


// labelHashSet patchSet;
// patchSet.insert(0);

// dist.internalField() = patchWave(mesh, patchSet, false).distance();


volScalarField AU
(
    IOobject
    (
        "AU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("AU", dimDensity/dimTime, 0)
);


label nWedgePatches = 0;
vector wedgeDirVec = vector::zero;
forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

        nWedgePatches++;
        wedgeDirVec += cmptMag(wpp.centreNormal());
    }
}

reduce(nWedgePatches, maxOp<label>());

if (nWedgePatches)
{
    Info<< nl << "Axisymmetric case: disabling the 3rd direction" << nl << endl;

    // We will const_cast as it is currently the tidiest way, until polyMesh is
    // modified or gives write access to solutionD
    Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

    reduce(wedgeDirVec, sumOp<vector>());

    wedgeDirVec /= mag(wedgeDirVec);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (wedgeDirVec[cmpt] > 1e-6)
        {
            solD[cmpt] = -1;
        }
        else
        {
            solD[cmpt] = 1;
        }
    }
}

