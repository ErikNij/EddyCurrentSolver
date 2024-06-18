
    vector totalSurfaceTensionForce =
        interface.totalSurfaceTensionForce();

    vector totalForce =
        interface.totalViscousForce()
      + interface.totalPressureForce();

    Info<< "Surface tension force: "
        << totalSurfaceTensionForce << endl;
    Info<< "Surface total force: "
        << totalForce << endl;

// TODO TEST: Sub-mesh
    const areaScalarField& K = interface.curvature();
//     const areaScalarField& K = interface.aMesh().faceCurvatures();

    Info<< "Surface curvature"
        << " min: " << gMin(K)
        << " max: " << gMax(K)
        << " average: " << gAverage(K)
        << endl;
