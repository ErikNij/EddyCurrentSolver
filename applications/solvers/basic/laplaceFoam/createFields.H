    Info<< "Reading field psi\n" << endl;

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Creating field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            word(),
            psi.dimensions()*dimLength,
            0
        ),
        calculatedFvPatchVectorField::typeName
    );


    Info<< "Creating field psiGrad\n" << endl;

    volVectorField psiGrad
    (
        IOobject
        (
            "psiGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(psi)
    );


    Info<< "Creating field psiGradFlux\n" << endl;

    volVectorField psiGradFlux
    (
        IOobject
        (
            "psiGradFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(psi)
    );


    Info<< "Creating field psiGradSnGrad\n" << endl;

    volVectorField psiGradSnGrad
    (
        IOobject
        (
            "psiGradSnGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(psi)
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity gamma\n" << endl;

    dimensionedScalar gamma
    (
        transportProperties.lookup("gamma")
    );

    dimensionedScalar source
    (
        transportProperties.lookup("source")
    );

    dimensionedScalar cellsource
    (
        source / sum(mesh.V()).value()
    );

    // Read solver for psi
    word solver
    (
        mesh.solutionDict().solver(psi.name()).lookup("solver")
    );

    bool deflated = false;
    if(solver == "deflatedPCG")
    {
        deflated = true;
    }

    label psiRefCell = 0;
    scalar psiRefValue = 0.0;
    setRefCell(psi, mesh.solutionDict().subDict("SIMPLE"), psiRefCell, psiRefValue);
    mesh.schemesDict().setFluxRequired(phi.name());
