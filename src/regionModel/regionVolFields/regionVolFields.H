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
    regionVolFields

Description
    ...

SourceFiles
    regionVolFields.C
    regionVolFieldMapping.C

\*---------------------------------------------------------------------------*/

#ifndef regionVolFields_H
#define regionVolFields_H

#include "fieldTypes.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "volMesh.H"

#include "regionGeometricField.H"
#include "regionGeoMesh.H"
#include "regionFvMesh.H"

// TODO [Low]: Implement template specializations for the currently
//             not-implemented membemer functions for boundary field mapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define declareRegionVolFieldMapping(Type)                                   \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvPatchField, volMesh, regionGeoMesh<regionFvMesh>                  \
>::mapInternalField                                                           \
(                                                                             \
    label regionI                                                             \
) const;                                                                      \
                                                                              \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvPatchField, volMesh, regionGeoMesh<regionFvMesh>                  \
>::rmapInternalField                                                          \
(                                                                             \
    label regionI                                                             \
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

declareRegionVolFieldMapping(scalar);
declareRegionVolFieldMapping(vector);
declareRegionVolFieldMapping(sphericalTensor);
declareRegionVolFieldMapping(symmTensor);
declareRegionVolFieldMapping(symmTensor4thOrder);
declareRegionVolFieldMapping(diagTensor);
declareRegionVolFieldMapping(tensor);


typedef regionGeometricField
    <scalar, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolScalarField;
typedef regionGeometricField
    <vector, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolVectorField;
typedef regionGeometricField
    <sphericalTensor, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolSphericalTensorField;
typedef regionGeometricField
    <symmTensor, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolSymmTensorField;
typedef regionGeometricField
    <symmTensor4thOrder, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolSymmTensor4thOrderField;
typedef regionGeometricField
    <diagTensor, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolDiagTensorField;
typedef regionGeometricField
    <tensor, fvPatchField, volMesh, regionGeoMesh<regionFvMesh> >
    regionVolTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
