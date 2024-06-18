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

#include "regionFvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector, Type>::type, fvPatchField, volMesh
//     >
// >
// grad
// (
//     const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
// )
// {
//     return fv::gaussGrad<Type>::gradf(ssf, "grad(" + ssf.name() + ')');
// }
//
//
// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector,Type>::type, fvPatchField, volMesh
//     >
// >
// grad
// (
//     const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
// )
// {
//     typedef typename outerProduct<vector, Type>::type GradType;
//     tmp<GeometricField<GradType, fvPatchField, volMesh> > Grad
//     (
//         fvc::grad(tssf())
//     );
//     tssf.clear();
//     return Grad;
// }


// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector,Type>::type, fvPatchField, volMesh
//     >
// >
// grad
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     const word& name
// )
// {
//     return fv::gradScheme<Type>::New
//     (
//         vf.mesh(),
//         vf.mesh().schemesDict().gradScheme(name)
//     )().grad(vf, name);
// }
//
//
// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector,Type>::type, fvPatchField, volMesh
//     >
// >
// grad
// (
//     const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
//     const word& name
// )
// {
//     tmp
//     <
//         GeometricField
//         <
//             typename outerProduct<vector, Type>::type, fvPatchField, volMesh
//         >
//     > tGrad
//     (
//         fvc::grad(tvf(), name)
//     );
//     tvf.clear();
//     return tGrad;
// }

template<class Type>
tmp
<
    regionGeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh,
        regionGeoMesh<regionFvMesh>
    >
> grad
(
    const regionGeometricField
    <
        Type, fvPatchField, volMesh,
        regionGeoMesh<regionFvMesh>
    >& vf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp
    <
        regionGeometricField
        <
            GradType, fvPatchField, volMesh,
            regionGeoMesh<regionFvMesh>
        >
    >
    tvfGrad
    (
        new regionGeometricField
        <
            GradType, fvPatchField, volMesh,
            regionGeoMesh<regionFvMesh>
        >
        (
            IOobject
            (
                "grad(" + vf.name() + ')',
                vf.mesh().time().timeName(),
                vf.mesh(),
                IOobject::NO_READ
            ),
            vf.mesh(),
            dimensioned<GradType>
            (
                word(),
                vf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    regionGeometricField
    <
        GradType, fvPatchField, volMesh,
        regionGeoMesh<regionFvMesh>
    >&
    vfgrad = tvfGrad();

    forAll(vf.mesh().regionNames(), regionI)
    {
        vfgrad[regionI] = fvc::grad(vf[regionI]);
    }

    return tvfGrad;
}

// template<class Type>
// tmp
// <
//     regionGeometricField
//     <
//         typename outerProduct<vector, Type>::type, fvPatchField, volMesh,
//         regionGeoMesh<regionFvMesh>
//     >
// > grad
// (
//     tmp
//     <
//         const regionGeometricField
//         <
//             Type, fvPatchField, volMesh,
//             regionGeoMesh<regionFvMesh>
//         >
//     >& tvf
// )
// {
//     return tvf;
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
