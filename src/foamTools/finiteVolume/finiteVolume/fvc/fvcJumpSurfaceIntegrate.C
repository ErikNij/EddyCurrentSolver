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

\*---------------------------------------------------------------------------*/

#include "fvcJumpSurfaceIntegrate.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    Field<Type>& vfIn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfNei
)
{
    const fvMesh& mesh = ssfOwn.mesh();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const Field<Type>& ssfOwnIn = ssfOwn.internalField();
    const Field<Type>& ssfNeiIn = ssfNei.internalField();

    forAll(owner, facei)
    {
        vfIn[owner[facei]] += ssfOwnIn[facei];
        vfIn[neighbour[facei]] -= ssfNeiIn[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const unallocLabelList& faceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& ssfOwnPatch =
            ssfOwn.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vfIn[faceCells[facei]] += ssfOwnPatch[facei];
        }
    }

    vfIn /= mesh.V();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfNei
)
{
    const fvMesh& mesh = ssfOwn.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssfOwn.name()+','+ssfNei.name()+')',
                ssfOwn.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssfOwn.dimensions()/dimVol,
                pTraits<Type>::zero
            ),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    surfaceIntegrate(vf.internalField(), ssfOwn, ssfNei);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssfOwn,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssfNei
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        fvc::surfaceIntegrate(tssfOwn(), tssfNei())
    );
    tssfOwn.clear();
    tssfNei.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssfNei
)
{
    const fvMesh& mesh = ssfOwn.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssfOwn.name()+','+ssfNei.name()+')',
                ssfOwn.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", ssfOwn.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();
    Field<Type>& vfIn = vf.internalField();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const Field<Type>& ssfOwnIn = ssfOwn.internalField();
    const Field<Type>& ssfNeiIn = ssfNei.internalField();

    forAll(owner, facei)
    {
        vfIn[owner[facei]] += ssfOwnIn[facei];
        vfIn[neighbour[facei]] += ssfNeiIn[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const unallocLabelList& faceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& ssfOwnPatch =
            ssfOwn.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            vfIn[faceCells[facei]] += ssfOwnPatch[facei];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssfOwn,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssfNei
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf =
        surfaceSum(tssfOwn(), tssfNei());
    tssfOwn.clear();
    tssfNei.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
