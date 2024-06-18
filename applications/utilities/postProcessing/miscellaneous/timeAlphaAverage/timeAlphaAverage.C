#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    Foam::timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);

    // Read region name
    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    // Create Mesh
    fvMesh mesh
    (
        IOobject
        (
            meshRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    volScalarField alpha1
    (
      IOobject
      (
        "alpha1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      ),
      mesh
    );

    volScalarField alpha1Mean
    (
      IOobject
      (
        "alpha1Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
      ),
      alpha1
    );

    volScalarField alpha1Diff
    (
      IOobject
      (
        "alpha1Diff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
      ),
      alpha1
    );


    alpha1Mean*=0;
    alpha1Diff*=0;

    int nfield=0;

    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (fieldHeader.headerOk())
        {

            mesh.readUpdate();

            volScalarField field(fieldHeader, mesh);
            alpha1Mean+=field;

            nfield++;

        }
        else
        {
            Info<< "    No field alpha1 " << endl;
        }

        Info<< endl;
    }

    alpha1Mean/=nfield;

    Info<< "writing alpha1Mean and asd to file at Time = "  << runTime.timeName() << endl << endl;

    alpha1Mean.write();

    forAll ( alpha1Diff,cellI )
    {
      alpha1Diff[cellI]=alpha1[cellI]-alpha1Mean[cellI];
    }

    alpha1Diff.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
