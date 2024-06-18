/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    decomposeRegionAddressing

Description
    Decompose region addressings after the case has been decomposed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "processorMeshes.H"

// TODO: Think about this agein. Currently we are "just" decomposing the region
//       addressings for all processors. But the labels are still addressing the
//       serial region! This is currently the INTENDED behaviour as it is needed
//       for the region-to-region addressing. Maybe it is not very intuitive and
//       should be revised together with the regionToRegionAddressing class!

// TODO: Warning/Exit if region option is not used (for defaultRegion). It
//       simply does not make any sense.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    // Determine the processor count directly
    label nProcs = 0;
    while (isDir(args.path()/(word("processor") + name(nProcs))))
    {
        ++nProcs;
    }

    if (!nProcs)
    {
        FatalErrorIn(args.executable())
            << "No processor* directories found"
            << exit(FatalError);
    }

    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );
    }

#   include "createNamedMesh.H"

    const word facesInstance = mesh.facesInstance();

    // Set all times on processor meshes equal to decomposed mesh
    forAll (databases, procI)
    {
        databases[procI].setTime(runTime);
    }

    // Read all meshes and addressing to reconstructed mesh
    processorMeshes procMeshes(databases, regionName);


    IOobject pointRegionAddObj
    (
        "pointRegionAddressing",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    IOobject pointMapObj
    (
        "pointMap",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!pointRegionAddObj.headerOk() && pointMapObj.headerOk())
    {
        pointRegionAddObj = pointMapObj;
    }


    if (pointRegionAddObj.headerOk())
    {
        Info<< "Decompose " << pointRegionAddObj.name() << " ... ";

        pointRegionAddObj.readOpt() = IOobject::MUST_READ;

        // Load addressing
        labelIOList regionPointAdd
        (
            pointRegionAddObj
        );

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const fvMesh& procMesh = procMeshes.meshes()[procI];

            const labelList& procPointAdd = procMeshes.pointProcAddressing()[procI];

            // Create proc region addressing
            labelIOList procRegionPointAdd
            (
                IOobject
                (
                    regionPointAdd.name(),
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                labelList(procPointAdd.size(), -1)
            );

            forAll (procRegionPointAdd, procPointI)
            {
                // Get point index of mesh
                label pointI = procPointAdd[procPointI];

                // Region addressing from region to default region
                procRegionPointAdd[procPointI] = regionPointAdd[pointI];
            }

            procRegionPointAdd.write();
        }

        Info<< "OK"
            << endl;
    }


    IOobject faceRegionAddObj
    (
        "faceRegionAddressing",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    IOobject faceMapObj
    (
        "faceMap",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!faceRegionAddObj.headerOk() && faceMapObj.headerOk())
    {
        faceRegionAddObj = faceMapObj;
    }

    if (faceRegionAddObj.headerOk())
    {
        Info<< "Decompose " << faceRegionAddObj.name() << " ... ";

        faceRegionAddObj.readOpt() = IOobject::MUST_READ;

        // Load addressing
        labelIOList regionFaceAdd
        (
            faceRegionAddObj
        );

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const fvMesh& procMesh = procMeshes.meshes()[procI];

            const labelList& procFaceAdd = procMeshes.faceProcAddressing()[procI];

            // Create proc region addressing
            labelIOList procRegionFaceAdd
            (
                IOobject
                (
                    regionFaceAdd.name(),
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                labelList(procFaceAdd.size(), -1)
            );

            forAll (procRegionFaceAdd, procFaceI)
            {
                // Regard face flip factors and shift to get face index
                // of mesh
                label faceI = mag(procFaceAdd[procFaceI]) - 1;

                // Keep face flip factors and shift in region addressing
                // from region to default region
                procRegionFaceAdd[procFaceI] = regionFaceAdd[faceI];
            }

            procRegionFaceAdd.write();
        }

        Info<< "OK"
            << endl;
    }


    IOobject cellRegionAddObj
    (
        "cellRegionAddressing",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    IOobject cellMapObj
    (
        "cellMap",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!cellRegionAddObj.headerOk() && cellMapObj.headerOk())
    {
        cellRegionAddObj = cellMapObj;
    }

    if (cellRegionAddObj.headerOk())
    {
        Info<< "Decompose " << cellRegionAddObj.name() << " ... ";

        cellRegionAddObj.readOpt() = IOobject::MUST_READ;

        // Load addressing
        labelIOList regionCellAdd
        (
            cellRegionAddObj
        );

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const fvMesh& procMesh = procMeshes.meshes()[procI];

            const labelList& procCellAdd = procMeshes.cellProcAddressing()[procI];

            // Create proc region addressing
            labelIOList procRegionCellAdd
            (
                IOobject
                (
                    regionCellAdd.name(),
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                labelList(procCellAdd.size(), -1)
            );

            forAll (procRegionCellAdd, procCellI)
            {
                // Get cell index of mesh
                label cellI = procCellAdd[procCellI];

                // Region addressing from region to default region
                procRegionCellAdd[procCellI] = regionCellAdd[cellI];
            }

            procRegionCellAdd.write();
        }

        Info<< "OK"
            << endl;
    }


    IOobject boundaryRegionAddObj
    (
        "boundaryRegionAddressing",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    IOobject patchMapObj
    (
        "patchMap",
        facesInstance,
        mesh.meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (!boundaryRegionAddObj.headerOk() && patchMapObj.headerOk())
    {
        boundaryRegionAddObj = patchMapObj;
    }

    if (boundaryRegionAddObj.headerOk())
    {
        Info<< "Decompose " << boundaryRegionAddObj.name() << " ... ";

        boundaryRegionAddObj.readOpt() = IOobject::MUST_READ;

        // Load addressing
        labelIOList regionBoundaryAdd
        (
            boundaryRegionAddObj
        );

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const fvMesh& procMesh = procMeshes.meshes()[procI];

            const labelList& procBoundaryAdd = procMeshes.boundaryProcAddressing()[procI];

            // Create proc region addressing
            labelIOList procRegionBoundaryAdd
            (
                IOobject
                (
                    regionBoundaryAdd.name(),
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                labelList(procBoundaryAdd.size(), -1)
            );

            forAll (procRegionBoundaryAdd, procBoundaryI)
            {
                // Get boundary index of mesh
                label boundaryI = procBoundaryAdd[procBoundaryI];

                // Skip non-existing boundaries of this proc
                if (boundaryI != -1)
                {
                    // Region addressing from region to default region
                    procRegionBoundaryAdd[procBoundaryI] = regionBoundaryAdd[boundaryI];
                }
            }

            procRegionBoundaryAdd.write();
        }

        Info<< "OK"
            << endl;
    }

    Info<< "End.\n" << endl;

    return(0);
}

// ************************************************************************* //
