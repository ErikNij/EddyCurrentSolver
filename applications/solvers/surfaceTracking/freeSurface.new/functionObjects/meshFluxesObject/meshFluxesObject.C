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

\*----------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>

#include "fvCFD.H"
#include "triPointRef.H"
#include "meshFluxesObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshFluxesObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        meshFluxesObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshFluxesObject::meshFluxesObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info << "Creating mesh fluxes check" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Triangle swept volume
scalar meshFluxesObject::triSweptVol
(
    const triPointRef& oT,
    const triPointRef& nT,
    const scalar deltaT,
    const label method
) const
{
    switch (method)
    {
        case 1:
        {
            scalar V1 =
            (
                ((nT.a() - oT.a()) & ((oT.b() - oT.a())^(oT.c() - oT.a())))
              + ((nT.b() - oT.b()) & ((oT.c() - oT.b())^(nT.a() - oT.b())))
              + ((oT.c() - nT.c()) & ((nT.b() - nT.c())^(nT.a() - nT.c())))
            );

            scalar V2 =
            (
                ((nT.a() - oT.a()) & ((oT.b() - oT.a())^(oT.c() - oT.a())))
              + ((oT.b() - nT.b()) & ((nT.a() - nT.b())^(nT.c() - nT.b())))
              + ((oT.c() - nT.c()) & ((oT.b() - nT.c())^(nT.a() - nT.c())))
            );

            // Current method (Frank Bos)
            return (1.0/12.0)*(V1 + V2);
        }

        case 2:
        {
            scalar V1 =
            (
                ((nT.a() - oT.a()) & ((oT.b() - oT.a())^(oT.c() - oT.a())))
              + ((nT.b() - oT.b()) & ((oT.c() - oT.b())^(nT.a() - oT.b())))
              + ((oT.c() - nT.c()) & ((nT.b() - nT.c())^(nT.a() - nT.c())))
            );

            scalar V2 =
            (
                ((nT.a() - oT.a()) & ((oT.b() - oT.a())^(oT.c() - oT.a())))
              + ((oT.b() - nT.b()) & ((nT.a() - nT.b())^(nT.c() - nT.b())))
              + ((oT.c() - nT.c()) & ((oT.b() - nT.c())^(nT.a() - nT.c())))
            );

            scalar V3 =
            (
                ((nT.b() - oT.b()) & ((oT.c() - oT.b())^(oT.a() - oT.b())))
              + ((oT.c() - nT.c()) & ((nT.b() - nT.c())^(nT.a() - nT.c())))
              + ((nT.a() - oT.a()) & ((nT.b() - oT.a())^(oT.c() - oT.a())))
            );

            scalar V4 =
            (
                ((nT.b() - oT.b()) & ((oT.c() - oT.b())^(oT.a() - oT.b())))
              + ((oT.a() - nT.a()) & ((nT.c() - nT.a())^(nT.b() - nT.a())))
              + ((oT.c() - nT.c()) & ((nT.b() - nT.c())^(oT.a() - nT.c())))
            );

            scalar V5 =
            (
                ((nT.c() - oT.c()) & ((oT.a() - oT.c())^(oT.b() - oT.c())))
              + ((oT.a() - nT.a()) & ((nT.c() - nT.a())^(nT.b() - nT.a())))
              + ((nT.b() - oT.b()) & ((nT.c() - oT.b())^(oT.a() - oT.b())))
            );

            scalar V6 =
            (
                ((nT.c() - oT.c()) & ((oT.a() - oT.c())^(oT.b() - oT.c())))
              + ((oT.b() - nT.b()) & ((nT.a() - nT.b())^(nT.c() - nT.b())))
              + ((oT.a() - nT.a()) & ((nT.c() - oT.a())^(oT.b() - nT.a())))
            );

            // Old method
            return
            (
                (1.0/6.0)*(1.0/6.0)*(V1 + V2 + V3 + V4 + V5 + V6)
            );
        }

        case 3:
        {
            // New Method (Blair Perot)
            return
            (
                ((nT.centre() - oT.centre())) &
                (
                    (0.5 * (nT.normal() + oT.normal()))
                  - (
                        (1.0 / 12.0) *
                        (
                            ((nT.a() - oT.a()) ^ (nT.b() - oT.b()))
                          + ((nT.b() - oT.b()) ^ (nT.c() - oT.c()))
                          + ((nT.c() - oT.c()) ^ (nT.a() - oT.a()))
                        )
                    )
                )
            );
        }

        default:
        {
            FatalErrorIn
            (
                "scalar meshFluxesObject::triSweptVol\n"
                "(\n"
                "    const triPointRef& oT,\n"
                "    const triPointRef& nT,\n"
                "    const label method\n"
                ") const\n"
            )
                << " Unknown method: " << method << nl
                << abort(FatalError);
        }
    }

    return 0.0;
}


//- Compute face swept-volume by central decomposition
scalar meshFluxesObject::sweptVol
(
    const face& f,
    const pointField& oldPoints,
    const pointField& newPoints,
    const scalar deltaT,
    const label method
) const
{
    if (method == 0)
    {
        return f.sweptVol(oldPoints, newPoints);
    }

    scalar sv = 0;

    if (f.size() == 3)
    {
        return triSweptVol
        (
            triPointRef
            (
                oldPoints[f[0]],
                oldPoints[f[1]],
                oldPoints[f[2]]
            ),
            triPointRef
            (
                newPoints[f[0]],
                newPoints[f[1]],
                newPoints[f[2]]
            ),
            deltaT,
            method
        );
    }

    point centreOldPoint = f.centre(oldPoints);
    point centreNewPoint = f.centre(newPoints);

    label nPoints = f.size();

    point nextOldPoint = centreOldPoint;
    point nextNewPoint = centreNewPoint;

    register label pI;

    for (pI = 0; pI < nPoints; pI++)
    {
        if (pI < nPoints - 1)
        {
            nextOldPoint = oldPoints[f[pI + 1]];
            nextNewPoint = newPoints[f[pI + 1]];
        }
        else
        {
            nextOldPoint = oldPoints[f[0]];
            nextNewPoint = newPoints[f[0]];
        }

        sv +=
        (
            triSweptVol
            (
                triPointRef
                (
                    centreOldPoint,
                    oldPoints[f[pI]],
                    nextOldPoint
                ),
                triPointRef
                (
                    centreNewPoint,
                    newPoints[f[pI]],
                    nextNewPoint
                ),
                deltaT,
                method
            )
        );
    }

    return sv;
}


//- Compute mesh-fluxes by method
tmp<scalarField> meshFluxesObject::sweptVols
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const pointField& oldPoints,
    const scalar deltaT,
    const label method
) const
{
    tmp<scalarField> tsweptVols(new scalarField(mesh.nFaces(), 0.0));

    scalarField& sweptVols = tsweptVols();

    // Create swept volumes
    const faceList& f = mesh.faces();

    forAll(f, faceI)
    {
        sweptVols[faceI] =
        (
            sweptVol
            (
                f[faceI],
                oldPoints,
                newPoints,
                deltaT,
                method
            )
        );
    }

    return tsweptVols;
}


// Set fluxes into GeometricField
void meshFluxesObject::setFluxes
(
    const scalar& rDeltaT,
    const scalarField& sweptVols,
    surfaceScalarField& phi
) const
{
    const fvMesh& mesh = phi.mesh();

    phi.internalField() =
    (
        scalarField::subField(sweptVols, mesh.nInternalFaces())
    );

    phi.internalField() *= rDeltaT;

    const fvPatchList& patches = mesh.boundary();

    forAll (patches, patchI)
    {
        phi.boundaryField()[patchI] = patches[patchI].patchSlice(sweptVols);
        phi.boundaryField()[patchI] *= rDeltaT;
    }
}


bool meshFluxesObject::start()
{
    return true;
}


bool meshFluxesObject::execute()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    Info << "Executing mesh fluxes function object." << endl;

    // Fetch time-step
    scalar deltaT = time_.deltaT().value();

    // Fetch old / new points
    const pointField& newPoints = mesh.points();
    const pointField& oldPoints = mesh.oldPoints();

    // Make the flux field
    surfaceScalarField meshPhi
    (
        IOobject
        (
            "MeshFluxes",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimVolume/dimTime, 0.0)
    );

    // Compute change in volume
    DimensionedField<scalar, volMesh> dV = (mesh.V() - mesh.V0());

    // Compute divMeshFlux by hand
    scalarField sv = sweptVols(mesh, newPoints, oldPoints, deltaT, 0)();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    scalarField divMeshPhi(mesh.nCells(), 0.0);

    forAll(own, faceI)
    {
        divMeshPhi[own[faceI]] += sv[faceI];
        divMeshPhi[nei[faceI]] -= sv[faceI];
    }

    forAll(mesh.boundary(), pi)
    {
        const unallocLabelList& pFaceCells = mesh.boundary()[pi].faceCells();

        label pStart = mesh.boundary()[pi].patch().start();

        forAll(pFaceCells, faceI)
        {
            divMeshPhi[pFaceCells[faceI]] += sv[pStart + faceI];
        }
    }

    Info << "Max error: " << max(dV - divMeshPhi) << endl;

    Info << "Mesh V: " << gSum(mesh.V())
         << " V0: " << gSum(mesh.V0())
         << endl;

    // The ddt term constructed by hand because it must be Euler
    DimensionedField<scalar, volMesh> dVdt =
    (
        (1.0 - (mesh.V0()/mesh.V())) / deltaT
    );

    // Update mesh-fluxes (with existing method in FOAM)
    setFluxes
    (
        (1.0 / deltaT),
        sweptVols(mesh, newPoints, oldPoints, deltaT, 0)(),
        meshPhi
    );

    volScalarField conserve = -fvc::div(meshPhi);
    conserve.internalField() += dVdt;

    Info << "Volume continuity errors: "
         << " max error = " << max(conserve.internalField())
         << endl;

    // Update mesh-fluxes (with current method)
    setFluxes
    (
        (1.0 / deltaT),
        sweptVols(mesh, newPoints, oldPoints, deltaT, 1)(),
        meshPhi
    );

    conserve = -fvc::div(meshPhi);
    conserve.internalField() += dVdt;

    Info << "Volume continuity errors (Bos): "
         << " max error = " << max(conserve.internalField())
         << endl;

    // Update mesh-fluxes (with old method)
    setFluxes
    (
        (1.0 / deltaT),
        sweptVols(mesh, newPoints, oldPoints, deltaT, 2)(),
        meshPhi
    );

    conserve = -fvc::div(meshPhi);
    conserve.internalField() += dVdt;

    Info << "Volume continuity errors (Jasak): "
         << " max error = " << max(conserve.internalField())
         << endl;

    // Update mesh-fluxes (with old method)
    setFluxes
    (
        (1.0 / deltaT),
        sweptVols(mesh, newPoints, oldPoints, deltaT, 3)(),
        meshPhi
    );

    conserve = -fvc::div(meshPhi);
    conserve.internalField() += dVdt;

    Info << "Volume continuity errors (Perot): "
         << " max error = " << max(conserve.internalField())
         << endl;

    return true;
}


bool meshFluxesObject::read(const dictionary& dict)
{
    return false;
}


} // End namespace Foam

// ************************************************************************* //
