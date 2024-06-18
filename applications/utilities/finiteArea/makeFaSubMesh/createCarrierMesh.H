
    fileName cMeshPrefix = "carrierMesh";
    if (mesh.name() != polyMesh::defaultRegion)
    {
        cMeshPrefix = mesh.name() + "_" + cMeshPrefix;
    }

    fvMeshSubset cMeshSubset
    (
        IOobject
        (
            cMeshPrefix,
            mesh.pointsInstance(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    cMeshSubset.setLargeCellSubset(aBMdata.faCells);
    fvMesh& cMesh = cMeshSubset.subMesh();
    Info << nl;
