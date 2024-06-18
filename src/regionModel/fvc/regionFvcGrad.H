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

InNamespace
    Foam::fvc

Description
    Calculate the gradient of the given field for region fields.

SourceFiles
    regionFvcGrad.C

\*---------------------------------------------------------------------------*/


#ifndef regionFvcGrad_H
#define regionFvcGrad_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

#include "fvcGrad.H"

#include "regionVolFields.H"

// TODO: Finish implementation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
//     template<class Type>
//     tmp
//     <
//         GeometricField
//         <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
//     > grad
//     (
//         const GeometricField<Type, fvsPatchField, surfaceMesh>&
//     );
//
//     template<class Type>
//     tmp
//     <
//         GeometricField
//         <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
//     > grad
//     (
//         const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
//     );

//     template<class Type>
//     tmp
//     <
//         GeometricField
//         <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
//     > grad
//     (
//         const GeometricField<Type, fvPatchField, volMesh>&,
//         const word& name
//     );
//
//     template<class Type>
//     tmp
//     <
//         GeometricField
//         <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
//     > grad
//     (
//         const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
//         const word& name
//     );

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
        >&
    );

//     template<class Type>
//     tmp
//     <
//         regionGeometricField
//         <
//             typename outerProduct<vector, Type>::type, fvPatchField, volMesh,
//             regionGeoMesh<regionFvMesh>
//         >
//     > grad
//     (
//         const tmp
//         <
//             regionGeometricField
//             <
//                 Type,
//                 fvPatchField,
//                 volMesh,
//                 regionGeoMesh<regionFvMesh>
//             >
//         >&
//     );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionFvcGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
