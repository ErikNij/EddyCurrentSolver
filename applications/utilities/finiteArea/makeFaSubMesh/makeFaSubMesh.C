/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    makeFaSubMesh

Description
    A mesh generator for a sub-polyMesh and corresponding sub-faMesh from a
    triangulation of a base-faMesh as part of base-polyMesh. Creation of
    the sub-polyMesh is based on a (1 cell) carrier mesh layer as fvMeshSubset
    of the base-polyMesh.

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "foamTime.H"
#include "argList.H"
#include "OSspecific.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "cellSet.H"
#include "triSurface.H"
#include "fvMeshSubset.H"
#include "Map.H"
#include "SortableList.H"
#include "Xfer.H"
#include "emptyPolyPatch.H"
#include "directMappedPolyPatch.H"
#include "directMappedWallPolyPatch.H"
#include "emptyFaPatch.H"
#include "demandDrivenData.H"

// TODO: Point ordering correct?

// TODO: Think about patchmap og sub-polyMesh: The new patch
//       'triangulatedFaces' comprises all faces from base-faMesh. So it may
//       consist of more then one patch from the base-polyMesh

// TODO: Point-, face-, edge- and patch-maps for sub-faMesh is missing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


struct meshData
{
    const fvMesh& mesh;
    int nInternalFaces;
    int nBoundaryFaces;
    int nActiveFaces;
    int nInactiveFaces;
    meshData
    (
        const fvMesh& mesh_
    )
    :
        mesh(mesh_),
        nInternalFaces(mesh.nInternalFaces()),
        nBoundaryFaces(-1),
        nActiveFaces(-1),
        nInactiveFaces(-1)
    {
        nBoundaryFaces = 0;
        forAll (mesh.boundaryMesh(), patchI)
        {
            nBoundaryFaces += mesh.boundaryMesh()[patchI].size();
        }

        nActiveFaces = nInternalFaces + nBoundaryFaces;

        nInactiveFaces = mesh.nFaces() - nActiveFaces;
    }
};


struct aMeshData
{
    const fvMesh& mesh;
    const faMesh& aMesh;
    const labelList& meshPoints;
    const labelList& faceLabels;
    labelList faceCells;
    cellSet faInternalCells;
    cellSet faBoundaryCells;
    cellSet faCells;
    aMeshData
    (
        const fvMesh& mesh_,
        const faMesh& aMesh_
    )
    :
        mesh(mesh_),
        aMesh(aMesh_),
        meshPoints(aMesh.patch().meshPoints()),
        faceLabels(aMesh.faceLabels()),
        faceCells(faceLabels.size(),-1),
        faInternalCells
        (
            mesh,
            "cellSet_faInternalCells",
            labelHashSet()
        ),
        faBoundaryCells
        (
            mesh,
            "cellSet_faBoundaryCells",
            labelHashSet()
        ),
        faCells
        (
            mesh,
            "cellSet_faCells",
            labelHashSet()
        )
    {
        // Calculate face cell addressing
        forAll (faceLabels, faceI)
        {
            label globalFaceI = faceLabels[faceI];

            label patchID = mesh.boundaryMesh().whichPatch(globalFaceI);
            label patchStart = mesh.boundaryMesh()[patchID].start();
            label patchFaceI = globalFaceI - patchStart;

            label patchFaceCellI =
                mesh.boundaryMesh()[patchID].faceCells()[patchFaceI];

            faceCells[faceI] = patchFaceCellI;
        }

        // Fill internal cell set
        faInternalCells.insert(faceCells);

        // Add all cells in contact with edges of the boundary
        // to the cell set faBoundaryCells
        forAll (aMesh.boundary(), patchI)
        {
            const faPatch& aPatch = aMesh.boundary()[patchI];

            // Get neighbour patch
            label ngbPolyPatchID = aPatch.ngbPolyPatchIndex();
            if (ngbPolyPatchID != -1)
            {
                labelList ngbPolyPatchFaces = aPatch.ngbPolyPatchFaces();

                forAll (aPatch, edgeI)
                {
                    label ngbGlobalFaceI = ngbPolyPatchFaces[edgeI];

                    label ngbPatchStart = mesh.boundaryMesh()[ngbPolyPatchID].start();
                    label ngbPatchFaceI = ngbGlobalFaceI - ngbPatchStart;

                    // ngbPatch cell
                    label ngpPatchFaceCellI =
                        mesh.boundaryMesh()[ngbPolyPatchID].faceCells()[ngbPatchFaceI];

                    faBoundaryCells.insert(ngpPatchFaceCellI);
                }
            }
        }

        // Remove all cells in contact with internal faces
        // from the cell set faBoundaryCells
        forAllConstIter (labelHashSet, faInternalCells, iter)
        {
            faBoundaryCells.unset(iter.key());
        }

        // Fill cell set
        faCells.insert(faInternalCells.toc());
        faCells.insert(faBoundaryCells.toc());

        // Write cellSets
        faInternalCells.Foam::regIOobject::write();
        faBoundaryCells.Foam::regIOobject::write();
        faCells.Foam::regIOobject::write();
    }
};


class aMeshTriangulation
{
private:
    const faMesh& aMesh_;
    pointField* points_;
    triFaceList* faces_;
    labelList* faceMap_;
    List<labelHashSet>* faceRmap_;
private:
    void calcPoints()
    {
        int nBasePoints = aMesh_.nPoints();
        int nCentrePoints = aMesh_.nFaces();
        int nPoints = nBasePoints + nCentrePoints;

        points_ = new pointField(nPoints, vector::zero);
        pointField& points = *points_;

        const pointField& basePoints = aMesh_.points();
        const pointField& centrePoints = aMesh_.patch().faceCentres();

        // Base vertices
        forAll(basePoints, basePointI)
        {
            label pointI = basePointI;

            points[pointI] = basePoints[pointI];
        }

        // Face centre vertices
        forAll(centrePoints, centrePointI)
        {
            points[nBasePoints + centrePointI] = centrePoints[centrePointI];
        }
    }

// TODO: Point ordering for points of triangulated faces
//       such that their normals point outward???
    void calcTriangulation()
    {

        int nBaseEdges = aMesh_.nEdges();
        int nBaseInternalEdges = aMesh_.nInternalEdges();
        int nFaces = nBaseEdges + nBaseInternalEdges;

        faces_ = new triFaceList(nFaces);
        triFaceList& faces = *faces_;

        faceMap_ = new labelList(nFaces, 0);
        labelList& faceMap = *faceMap_;

        faceRmap_ = new List<labelHashSet>(aMesh_.nFaces(), labelHashSet());
        List<labelHashSet>& faceRmap = *faceRmap_;

        const edgeList& baseEdges = aMesh_.edges();

        // Internal triangulation based on lower/upper addressing
        int nBasePoints = aMesh_.nPoints();
        const unallocLabelList& baseOwn = aMesh_.owner();
        const unallocLabelList& baseNei = aMesh_.neighbour();

        forAll(baseOwn, baseEdgeI)
        {
            edge e = baseEdges[baseEdgeI];

            // Owner side triangle
            label ownP1 = e.start();
            label ownP2 = e.end();
            label ownP3 = nBasePoints + baseOwn[baseEdgeI];

            label ownFaceI = 2*baseEdgeI;

            faces[ownFaceI] = triFace(ownP1, ownP2, ownP3);
            faceMap[ownFaceI] = baseOwn[baseEdgeI];
            faceRmap[baseOwn[baseEdgeI]].insert(ownFaceI);

            // Neighbour side triangle
            label neiP1 = e.end();
            label neiP2 = e.start();
            label neiP3 = nBasePoints + baseNei[baseEdgeI];

            label neiFaceI = 2*baseEdgeI + 1;

            faces[neiFaceI] = triFace(neiP1, neiP2, neiP3);
            faceMap[neiFaceI] = baseNei[baseEdgeI];
            faceRmap[baseNei[baseEdgeI]].insert(neiFaceI);
        }

        // Patch boundary triangulation
        const faPatchList& basePatches = aMesh_.boundary();
        int nBasePatchEdges = 0;
        forAll(basePatches, basePatchI)
        {
            const faPatch& basePatch = basePatches[basePatchI];

            if (basePatch.size() > 0)
            {
                const unallocLabelList& basePatchEdgeFaces =
                    basePatch.edgeFaces();

                forAll(basePatch, basePatchEdgeI)
                {
                    // Start corresponding to current patch
                    label baseEdgeI = basePatch.start() + basePatchEdgeI;

                    edge e = baseEdges[baseEdgeI];

                    // Patch edge triangle
                    label pP1 = e.start();
                    label pP2 = e.end();
                    label pP3 =
                        nBasePoints + basePatchEdgeFaces[basePatchEdgeI];

                    label pFaceI =
                        2*nBaseInternalEdges + nBasePatchEdges++;

                    faces[pFaceI] = triFace(pP1, pP2, pP3);
                    faceMap[pFaceI] = basePatchEdgeFaces[basePatchEdgeI];
                    faceRmap[basePatchEdgeFaces[basePatchEdgeI]].insert(pFaceI);
                }
            }
        }

        // Empty boundary triangulation
        int nBaseNonEmptyEdges = nBaseInternalEdges + nBasePatchEdges;
        int nBaseEmptyEdges = nBaseEdges - nBaseNonEmptyEdges;

        if (nBaseEmptyEdges > 0)
        {
            labelList::subList baseEmptyOwn
            (
                aMesh_.edgeOwner(),
                nBaseEmptyEdges,
                nBaseNonEmptyEdges
            );

            forAll(baseEmptyOwn, baseEmptyEdgeI)
            {
                // Start after all non-empty edges
                label baseEdgeI = nBaseNonEmptyEdges + baseEmptyEdgeI;

                edge e = baseEdges[baseEdgeI];

                // Empty edge triangle
                label eP1 = e.start();
                label eP2 = e.end();
                label eP3 = nBasePoints + baseEmptyOwn[baseEmptyEdgeI];

                label eFaceI =
                    2*nBaseInternalEdges + nBasePatchEdges + baseEmptyEdgeI;

                faces[eFaceI] = triFace(eP1, eP2, eP3);
                faceMap[eFaceI] = baseEmptyOwn[baseEmptyEdgeI];
                faceRmap[baseEmptyOwn[baseEmptyEdgeI]].insert(eFaceI);
            }
        }
    }
public:
    aMeshTriangulation
    (
        const faMesh& aMesh
    )
    :
        aMesh_(aMesh),
        points_(NULL),
        faces_(NULL),
        faceMap_(NULL),
        faceRmap_(NULL)
    {
        calcPoints();
        calcTriangulation();
    }

    ~aMeshTriangulation()
    {
        deleteDemandDrivenData(faceRmap_);
        deleteDemandDrivenData(faceMap_);
        deleteDemandDrivenData(faces_);
        deleteDemandDrivenData(points_);
    }

    const pointField& points() const
    {
        return *points_;
    }

    const triFaceList& faces() const
    {
        return *faces_;
    }

    const labelList& faceMap() const
    {
        return *faceMap_;
    }

    const List<labelHashSet>& faceRmap() const
    {
        return *faceRmap_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"
#   include "createFaMesh.H"

    Switch debug(false);

    // Base mesh data
    meshData bMdata(mesh);

    if (debug)
    {
        Info << "DEBUG | bMdata.nInternalFaces = " << bMdata.nInternalFaces << endl;
        Info << "DEBUG | bMdata.nBoundaryFaces = " << bMdata.nBoundaryFaces << endl;
        Info << "DEBUG | bMdata.nActiveFaces = " << bMdata.nActiveFaces << endl;
        Info << "DEBUG | bMdata.nInactiveFaces = " << bMdata.nInactiveFaces << endl;
        Info << endl;
    }

    // Area base-mesh data
    aMeshData aBMdata(mesh, aMesh);

    if (debug)
    {
        Info << "DEBUG | aBMdata.internalCells.size() = " << aBMdata.faInternalCells.size() << endl;
        Info << "DEBUG | aBMdata.boundaryCells.size() = " << aBMdata.faBoundaryCells.size() << endl;
        Info << "DEBUG | aBMdata.cells.size() = " << aBMdata.faCells.size() << endl;
        Info << endl;
    }

    // Create carrier-mesh
#   include "createCarrierMesh.H"

    // Carrier-mesh data
    meshData cMdata(cMesh);

    if (debug)
    {
        Info << "DEBUG | cMdata.nInternalFaces = " << cMdata.nInternalFaces << endl;
        Info << "DEBUG | cMdata.nBoundaryFaces = " << cMdata.nBoundaryFaces << endl;
        Info << "DEBUG | cMdata.nActiveFaces = " << cMdata.nActiveFaces << endl;
        Info << "DEBUG | cMdata.nInactiveFaces = " << cMdata.nInactiveFaces << endl;
        Info << endl;
    }

    // Area base-mesh triangulation
    aMeshTriangulation aBMtri(aMesh);

    if (debug)
    {
        Info << "DEBUG | aBMtri.points().size() = " << aBMtri.points().size() << endl;
        Info << "DEBUG | aBMtri.faces().size() = " << aBMtri.faces().size() << endl;
        Info << endl;
    }

// Points

    // Point map from base-mesh to carrier-mesh
    // Key: Point label in base-mesh
    // Content: Point label in carrier-mesh
    Map<label> bMpointsCMmap;
    Map<label> usedBMpointsCMmap;

    // Collect all base-mesh point labels from points in carrier-mesh
    forAll(cMeshSubset.pointMap(), pointI)
    {
            bMpointsCMmap.insert(cMeshSubset.pointMap()[pointI], pointI);
    }
    usedBMpointsCMmap = bMpointsCMmap;

    // Remove all base-mesh point labels of area base-mesh
    forAll(aBMdata.meshPoints, pointI)
    {
            usedBMpointsCMmap.erase(aBMdata.meshPoints[pointI]);
    }

    // Point map from carrier-mesh to base-mesh
    // Key: Point label in carrier-mesh
    // Content: Point label in base-mesh
    Map<label> usedCMpointsBMmap;

    forAllConstIter (Map<label>, usedBMpointsCMmap, iter)
    {
        label cMpointI = iter();

        usedCMpointsBMmap.insert(cMpointI, cMeshSubset.pointMap()[cMpointI]);
    }

    // Init points
    int nPoints = usedCMpointsBMmap.size() + aBMtri.points().size();
    pointField points(nPoints, vector::zero);
    labelList usedCMpoints = usedCMpointsBMmap.sortedToc();

    // Add used points of carrier-mesh
    forAll(usedCMpoints, pointI)
    {
        label cMpointI = usedCMpoints[pointI];

        points[pointI] = cMesh.points()[cMpointI];
    }

    // Add all points of triangulated area base-mesh
    forAll (aBMtri.points(), pointI)
    {
        points[usedCMpoints.size() + pointI] = aBMtri.points()[pointI];
    }

    // Create reverse pointMap for carrier-mesh
    Map<label> cMpointRmap;
    forAll (usedCMpoints, pointI)
    {
        label cMpointI = usedCMpoints[pointI];

        cMpointRmap.insert
        (
            cMpointI,
            pointI
        );
    }

    // Create reverse pointMap for area base-mesh
    Map<label> aBMpointRmap;
    forAll (aBMdata.meshPoints, pointI)
    {
        aBMpointRmap.insert
        (
            bMpointsCMmap[aBMdata.meshPoints[pointI]],
            usedCMpoints.size() + pointI
        );
    }

    if (debug)
    {
        Info << "DEBUG | usedCMpoints.size() = " << usedCMpoints.size() << endl;
        Info << "DEBUG | points.size() = " << points.size() << endl;
        Info << endl;
    }

// Faces

    // Face map from base-mesh to carrier-mesh
    // Key: Face label in base-mesh
    // Content: Face label in carrier-mesh
    Map<label> bMfacesCMmap;
    Map<label> usedBMfacesCMmap;

    // Collect all base-mesh face labels from faces in carrier-mesh
    forAll(cMeshSubset.faceMap(), faceI)
    {
            bMfacesCMmap.insert(cMeshSubset.faceMap()[faceI], faceI);
    }
    usedBMfacesCMmap = bMfacesCMmap;

    // Remove all base-mesh face labels of area base-mesh
    forAll(aBMdata.faceLabels, faceI)
    {
            usedBMfacesCMmap.erase(aBMdata.faceLabels[faceI]);
    }

    // Face map from carrier-mesh to base-mesh
    // Key: Face label in carrier-mesh
    // Content: Face label in base-mesh
    Map<label> usedCMfacesBMmap;

    forAllConstIter (Map<label>, usedBMfacesCMmap, iter)
    {
        label faceI = iter();

        usedCMfacesBMmap.insert(faceI, cMeshSubset.faceMap()[faceI]);
    }

    // Init faces
    int nFaces = usedCMfacesBMmap.size() + aBMtri.faces().size();
    faceList faces(nFaces);
    labelList usedCMfaces = usedCMfacesBMmap.sortedToc();

    // Add used faces of carrier-mesh
    forAll(usedCMfaces, faceI)
    {
        label cMfaceI = usedCMfaces[faceI];

        face curFace = cMesh.faces()[cMfaceI];

        // Map face vertices
        labelList pointLabels(curFace.size(), -1);
        forAll(pointLabels, pointI)
        {
            // Point label is part of usedBMpoints (point from carrier-mesh)
            if (usedCMpointsBMmap.found(curFace[pointI]))
            {
                pointLabels[pointI] = cMpointRmap[curFace[pointI]];
            }
            // Point label is NOT part of usedBMpoints (point from area base-mesh)
            else
            {
                pointLabels[pointI] = aBMpointRmap[curFace[pointI]];
            }
        }

        faces[faceI] = face(pointLabels);
    }

    // Add all faces of triangulated area base-mesh
    forAll (aBMtri.faces(), faceI)
    {
        face curFace = aBMtri.faces()[faceI];

        // Map face vertices
        labelList pointLabels(curFace.size(), -1);
        forAll(pointLabels, pointI)
        {
            pointLabels[pointI] = usedCMpoints.size() + curFace[pointI];
        }

        faces[usedCMfaces.size() + faceI] = face(pointLabels);
    }

    if (debug)
    {
        Info << "DEBUG | usedCMfaces.size() = " << usedCMfaces.size() << endl;
        Info << "DEBUG | faces.size() = " << faces.size() << endl;
        Info << endl;
    }

// Patches

    // Init patch face data
    label nInternalFaces = 0;
    label nBoundaryFaces = 0;

    // Init patch data
    label patchI = -1;
    label lastPatchID = -1;
    label lastPatchStart = 0;
    Map<label> patchLabels;
    Map<word> patchTypes;
    Map<word> patchNames;
    Map<label> patchSizes;
    Map<label> patchStarts;

    // Count patches and gather face information for carrier-mesh
    forAll (usedCMfaces, faceI)
    {
        label cMfaceI = usedCMfaces[faceI];

        label curPatchID = cMesh.boundaryMesh().whichPatch(cMfaceI);

        if (curPatchID < 0)
        {
            nInternalFaces++;
        }
        else
        {
            nBoundaryFaces++;

            if (lastPatchID != curPatchID)
            {
                patchI++;
                patchLabels.insert(patchI, curPatchID);
                patchTypes.insert(patchI, cMesh.boundary()[curPatchID].type());
                patchNames.insert(patchI, cMesh.boundary()[curPatchID].name());
                if (patchI > 0)
                {
                    patchSizes.insert(patchI-1, faceI - lastPatchStart);
                }
                patchStarts.insert(patchI, faceI);

                lastPatchID = curPatchID;
                lastPatchStart = faceI;
            }
        }
    }
    patchSizes.insert(patchI, usedCMfaces.size() - lastPatchStart);

    // Add/Modify patch data for triangulated area base-mesh
    nBoundaryFaces += aBMtri.faces().size();
    ++patchI;
// TODO: Which patch label(s) for triangulatedFaces? May be more then one patch...
    patchLabels.insert(patchI, -1);
    patchTypes.insert(patchI, "patch");
    patchNames.insert(patchI, "triangulatedFaces");
    patchStarts.insert(patchI, usedCMfaces.size());
    patchSizes.insert(patchI, aBMtri.faces().size());

    if (debug)
    {
        Info << "DEBUG | nInternalFaces = " << nInternalFaces << endl;
        Info << "DEBUG | nBoundaryFaces = " << nBoundaryFaces << endl;
        Info << "DEBUG | patchTypes = " << patchTypes << endl;
        Info << "DEBUG | patchNames = " << patchNames << endl;
        Info << "DEBUG | patchSizes = " << patchSizes << endl;
        Info << "DEBUG | patchStarts = " << patchStarts << endl;
        Info << endl;
    }

// Owner/Neighbour

    // Face owner
    labelList faceOwner(nFaces, -1);
    forAll(usedCMfaces, faceI)
    {
        label cMfaceI = usedCMfaces[faceI];

        faceOwner[faceI] = cMesh.faceOwner()[cMfaceI];
    }
    forAll (aBMtri.faces(), faceI)
    {
        label aBMfaceI = aBMtri.faceMap()[faceI];
        label bMfaceI = aBMdata.faceLabels[aBMfaceI];
        label cMfaceI = bMfacesCMmap[bMfaceI];

        faceOwner[usedCMfaces.size() + faceI] = cMesh.faceOwner()[cMfaceI];
    }

    // Face neighbour
    labelList faceNeighbour(nInternalFaces, -1);
    forAll(faceNeighbour, faceI)
    {
        label cMfaceI = usedCMfaces[faceI];

        faceNeighbour[faceI] = cMesh.faceNeighbour()[cMfaceI];
    }

    if (debug)
    {
        Info << "DEBUG | faceOwner.size() = " << faceOwner.size() << endl;
        Info << "DEBUG | faceNeighbour.size() = " << faceNeighbour.size() << endl;
        Info << endl;
    }


// subPolyMesh

    // Creation

        polyMesh sMesh
        (
            IOobject
            (
                mesh.name() + "_faSubMesh",
                mesh.pointsInstance(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(points),
            xferCopy(faces),
            xferCopy(faceOwner),
            xferCopy(faceNeighbour),
            false
        );

        // Create patches
        List<polyPatch*> sPatches(patchTypes.size(), NULL);
        forAll (sPatches, patchI)
        {
            word curPatchType = patchTypes[patchI];

            // Empty directions may become invalidated as we have
            // triangulated faces
            if (curPatchType == emptyPolyPatch::typeName)
            {
                curPatchType = polyPatch::typeName;
            }

            // Mapped region will not be present
            if (curPatchType == directMappedPolyPatch::typeName)
            {
                curPatchType = polyPatch::typeName;
            }

            // Mapped region will not be present
            if (curPatchType == directMappedWallPolyPatch::typeName)
            {
                curPatchType = wallPolyPatch::typeName;
            }

            sPatches[patchI] = polyPatch::New
            (
                curPatchType,
                patchNames[patchI],
                patchSizes[patchI],
                patchStarts[patchI],
                patchI,
                sMesh.boundaryMesh()
            ).ptr();
        }

        // Add patches
        sMesh.addPatches(sPatches, true);

        // Write subMesh
        Info << "Write subMesh ... ";
        sMesh.write();
        Info << "Done" << endl << endl;


    // Maps

        {
            labelIOList pointMap
            (
                IOobject
                (
                    "pointMap",
                    sMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    sMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                labelList(points.size(), -1)
            );
            forAll (usedCMpoints, pointI)
            {
                label cMpointI = usedCMpoints[pointI];
                label bMpointI = usedCMpointsBMmap[cMpointI];

                pointMap[pointI] = bMpointI;
            }
            forAll (aBMdata.meshPoints, pointI)
            {
                pointMap[usedCMpoints.size() + pointI] = aBMdata.meshPoints[pointI];
            }
            forAll (aBMdata.faceLabels, pointI)
            {
                pointMap[usedCMpoints.size() + aBMdata.meshPoints.size() + pointI] = -1;
            }

            if (mesh.name() != polyMesh::defaultRegion)
            {
                IOobject pointRegionAddressingObj
                (
                    "pointRegionAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                IOobject pointRegionMapObj
                (
                    "pointMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                if (!pointRegionAddressingObj.headerOk() && pointRegionMapObj.headerOk())
                {
                    pointRegionAddressingObj = pointRegionMapObj;
                }

                pointRegionAddressingObj.readOpt() = IOobject::MUST_READ;

                labelIOList pointRegionAdd = labelIOList(pointRegionAddressingObj);

                forAll (pointMap, pointI)
                {
                    if (pointMap[pointI] != -1)
                    {
                        pointMap[pointI] = pointRegionAdd[pointMap[pointI]];
                    }
                }
            }

            Info << "Write pointMap ... ";
            pointMap.write();
            Info << "Done" << endl;
        }

        {
            labelIOList faceMap
            (
                IOobject
                (
                    "faceMap",
                    sMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    sMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                labelList(faces.size(), -1)
            );
            forAll (usedCMfaces, faceI)
            {
                label cMfaceI = usedCMfaces[faceI];
                label bMfaceI = usedCMfacesBMmap[cMfaceI];

                faceMap[faceI] = bMfaceI;
            }
            forAll (aBMtri.faces(), faceI)
            {
                faceMap[usedCMfaces.size() + faceI] = aBMdata.faceLabels[aBMtri.faceMap()[faceI]];
            }

            if (mesh.name() != polyMesh::defaultRegion)
            {
                IOobject faceRegionAddressingObj
                (
                    "faceRegionAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                IOobject faceRegionMapObj
                (
                    "faceMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                bool fromMap = false;

                if (!faceRegionAddressingObj.headerOk() && faceRegionMapObj.headerOk())
                {
                    faceRegionAddressingObj = faceRegionMapObj;

                    fromMap = true;
                }

                faceRegionAddressingObj.readOpt() = IOobject::MUST_READ;

                labelIOList faceRegionAdd = labelIOList(faceRegionAddressingObj);

                if (!fromMap)
                {
                    forAll (faceRegionAdd, faceI)
                    {
                        label magFaceRegionMap = mag(faceRegionAdd[faceI]);

                        faceRegionAdd[faceI] = magFaceRegionMap - 1;
                    }
                }

                forAll (faceMap, faceI)
                {
                    if (faceMap[faceI] != -1)
                    {
                        faceMap[faceI] = faceRegionAdd[faceMap[faceI]];
                    }
                }
            }

            Info << "Write faceMap ... ";
            faceMap.write();
            Info << "Done" << endl;
        }

        {
            labelIOList cellMap
            (
                IOobject
                (
                    "cellMap",
                    sMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    sMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                cMeshSubset.cellMap()
            );

            if (mesh.name() != polyMesh::defaultRegion)
            {
                IOobject cellRegionAddressingObj
                (
                    "cellRegionAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                IOobject cellRegionMapObj
                (
                    "cellMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                if (!cellRegionAddressingObj.headerOk() && cellRegionMapObj.headerOk())
                {
                    cellRegionAddressingObj = cellRegionMapObj;
                }

                cellRegionAddressingObj.readOpt() = IOobject::MUST_READ;

                labelIOList cellRegionAdd = labelIOList(cellRegionAddressingObj);

                forAll (cellMap, cellI)
                {
                    if (cellMap[cellI] != -1)
                    {
                        cellMap[cellI] = cellRegionAdd[cellMap[cellI]];
                    }
                }
            }

            Info << "Write cellMap ... ";
            cellMap.write();
            Info << "Done" << endl;
        }

        {
            labelIOList patchMap
            (
                IOobject
                (
                    "patchMap",
                    sMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    sMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                labelList(patchLabels.size(), -1)
            );
            forAll (patchMap, patchI)
            {
                patchMap[patchI] = patchLabels[patchI];
            }

            if (mesh.name() != polyMesh::defaultRegion)
            {
                IOobject boundaryRegionAddressingObj
                (
                    "boundaryRegionAddressing",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                IOobject patchRegionMapObj
                (
                    "patchMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                );

                if (!boundaryRegionAddressingObj.headerOk() && patchRegionMapObj.headerOk())
                {
                    boundaryRegionAddressingObj = patchRegionMapObj;
                }

                boundaryRegionAddressingObj.readOpt() = IOobject::MUST_READ;

                labelIOList boundaryRegionAdd = labelIOList(boundaryRegionAddressingObj);

                forAll (patchMap, patchI)
                {
                    if (patchMap[patchI] != -1)
                    {
                        patchMap[patchI] = boundaryRegionAdd[patchMap[patchI]];
                    }
                }
            }

            Info << "Write patchMap ... ";
            patchMap.write();
            Info << "Done" << endl;
        }


// subAreaMesh

    // faMeshDefinition

        // Read faMeshDefinition dictionary from area base-mesh
        IOdictionary aBMdefinition
        (
            IOobject
            (
                "faMeshDefinition",
                mesh.facesInstance(),
                faMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Create faMeshDefinition dictionary
        IOdictionary faMeshDefinition
        (
            IOobject
            (
                "faMeshDefinition",
                sMesh.facesInstance(),
                faMesh::meshSubDir,
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        // Set polyMeshPatches
        faMeshDefinition.add
        (
            "polyMeshPatches",
            wordList(1, "triangulatedFaces")
        );

        // Copy boundary from area base-mesh
        faMeshDefinition.add
        (
            "boundary",
            aBMdefinition.subDict("boundary")
        );

        dictionary& bndDict = faMeshDefinition.subDict("boundary");

        wordList faPatchNames = bndDict.toc();

        // Check for empty patches and replace all
        // ownerPolyPatches with triangulatedFaces
        forAll (faPatchNames, patchI)
        {
            dictionary& curPatchDict =
                bndDict.subDict(faPatchNames[patchI]);

            word curPatchType = word(curPatchDict.lookup("type"));

            // Empty directions may become invalidated as we have
            // triangulated faces
            if (curPatchType == emptyFaPatch::typeName)
            {
                curPatchType = faPatch::typeName;
            }

            curPatchDict.set("type", curPatchType);

            curPatchDict.set("ownerPolyPatch", word("triangulatedFaces"));
        }

        // Write faMeshDefinition
        faMeshDefinition.regIOobject::write();

    // Creation

        faMesh sAmesh(sMesh, fileName("faMeshDefinition"));
        sAmesh.write();

    // Maps

// TODO: Point-, face-, edge- and patch-maps

    Info<< endl << "End.\n" << endl;

    return (0);
}

// ************************************************************************* //

