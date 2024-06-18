
    // Update gradient of V in default region
    // from data in conductor region
    if (control.meshIs3D())
    {
        if (Control::debug)
        {
            Info<< Control::typeName << " | VGradUpdate.H : "
                << "Update gradient of V."
                << endl;
        }

// TEST: Use relaxation?
//         globalStorage.VReGrad()[DEFAULT].storePrevIter();
//         globalStorage.VImGrad()[DEFAULT].storePrevIter();

        // Update gradient of V
        globalStorage.VReGrad().rmap(CONDUCTOR);
        globalStorage.VImGrad().rmap(CONDUCTOR);

// TEST: Use relaxation?
//         const scalar VReGradRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.VReGrad()[DEFAULT].name());
//         const scalar VImGradRf = globalMesh[DEFAULT].solutionDict().fieldRelaxationFactor(globalStorage.VImGrad()[DEFAULT].name());
//         globalStorage.VReGrad()[DEFAULT].relax(VReGradRf);
//         globalStorage.VImGrad()[DEFAULT].relax(VImGradRf);
    }
