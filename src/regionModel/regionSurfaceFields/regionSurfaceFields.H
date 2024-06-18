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

Class
    regionSurfaceFields

Description
    ...

SourceFiles
    regionSurfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef regionSurfaceFields_H
#define regionSurfaceFields_H

#include "fieldTypes.H"
#include "surfaceFields.H"
#include "fvsPatchFields.H"
#include "surfaceMesh.H"

#include "regionGeometricField.H"
#include "regionGeoMesh.H"
#include "regionFvMesh.H"

// TODO [Low]: Implement template specializations for the currently
//             not-implemented membemer functions for boundary field mapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define declareRegionSurfaceFieldMapping(Type)                               \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh>             \
>::mapInternalField                                                           \
(                                                                             \
    label regionI                                                             \
) const;                                                                      \
                                                                              \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh>             \
>::rmapInternalField                                                          \
(                                                                             \
    label regionI                                                             \
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef regionGeometricField
    <scalar, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh> >
    regionSurfaceScalarField;
typedef regionGeometricField
    <vector, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh> >
    regionSurfaceVectorField;
typedef regionGeometricField
    <tensor, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh> >
    regionSurfaceTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
