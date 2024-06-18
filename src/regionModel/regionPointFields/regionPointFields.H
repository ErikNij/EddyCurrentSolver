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
    regionPointFields

Description
    ...

SourceFiles
    regionPointFields.C

\*---------------------------------------------------------------------------*/

#ifndef regionPointFields_H
#define regionPointFields_H

#include "fieldTypes.H"
#include "pointFields.H"
#include "pointPatchFields.H"
#include "pointMesh.H"

#include "regionGeometricField.H"
#include "regionGeoMesh.H"
#include "regionPointMesh.H"

// TODO [Low]: Implement template specializations for the currently
//             not-implemented membemer functions for boundary field mapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define declareRegionSurfaceFieldMapping(Type)                               \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, pointPatchField, pointMesh, regionGeoMesh<regionPointMesh>          \
>::mapInternalField                                                           \
(                                                                             \
    label regionI                                                             \
) const;                                                                      \
                                                                              \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, pointPatchField, pointMesh, regionGeoMesh<regionPointMesh>          \
>::rmapInternalField                                                          \
(                                                                             \
    label regionI                                                             \
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef regionGeometricField
    <scalar, pointPatchField, pointMesh, regionGeoMesh<regionPointMesh> >
    regionPointScalarField;
typedef regionGeometricField
    <vector, pointPatchField, pointMesh, regionGeoMesh<regionPointMesh> >
    regionPointVectorField;
typedef regionGeometricField
    <tensor, pointPatchField, pointMesh, regionGeoMesh<regionPointMesh> >
    regionPointTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
