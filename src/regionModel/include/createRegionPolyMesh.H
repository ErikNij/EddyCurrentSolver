//
// createRegionPolyMesh.H
// ~~~~~~~~~~~~~~~~~~

    Foam::Info
        << "Create region polyMesh for time = "
        << runTime.timeName() << endl;

    Foam::regionPolyMesh regionMesh(runTime);

    Foam::Info << Foam::endl;
