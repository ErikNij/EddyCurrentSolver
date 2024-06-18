//
// createRegionMeshUninitialized.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionFvMesh regionMesh(runTime, false);

    Foam::Info << Foam::endl;
