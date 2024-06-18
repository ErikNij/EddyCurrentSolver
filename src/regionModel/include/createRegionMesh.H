//
// createRegionMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionFvMesh regionMesh(runTime);

    Foam::Info << Foam::endl;
