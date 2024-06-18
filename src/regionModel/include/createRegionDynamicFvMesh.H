//
// createRegionDynamicFvMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create dynamic region mesh for time = "
        << runTime.timeName() << endl;

    Foam::regionDynamicFvMesh regionMesh(runTime);

    Foam::Info << Foam::endl;
