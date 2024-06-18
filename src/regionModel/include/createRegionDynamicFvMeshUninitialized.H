//
// createRegionDynamicFvMeshUninitialized.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create dynamic region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionDynamicFvMesh regionMesh(runTime, false);

    Foam::Info << Foam::endl;
