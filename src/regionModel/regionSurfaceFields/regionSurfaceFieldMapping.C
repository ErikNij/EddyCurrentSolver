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

#include "regionSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define defineRegionSurfaceFieldMapping(Type)                                \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh>             \
>::mapInternalField                                                           \
(                                                                             \
    label regionI                                                             \
) const                                                                       \
{                                                                             \
    label regionI0 = regions()[polyMesh::defaultRegion];                      \
                                                                              \
    const GeometricFieldType& vf0 = operator[](regionI0);                     \
                                                                              \
    GeometricFieldType& vf = operator[](regionI);                             \
                                                                              \
    const labelList& map =                                                    \
        mesh().typeMap(addressingTypes::FACE, regionI, regionI0);             \
                                                                              \
    forAll (vf, facei)                                                        \
    {                                                                         \
        vf[facei] = vf0[map[facei]];                                          \
    }                                                                         \
}                                                                             \
                                                                              \
                                                                              \
template<>                                                                    \
void regionGeometricField                                                     \
<                                                                             \
    Type, fvsPatchField, surfaceMesh, regionGeoMesh<regionFvMesh>             \
>::rmapInternalField                                                          \
(                                                                             \
    label regionI                                                             \
) const                                                                       \
{                                                                             \
    label regionI0 = regions()[polyMesh::defaultRegion];                      \
                                                                              \
    GeometricFieldType& vf0 = operator[](regionI0);                           \
                                                                              \
    const GeometricFieldType& vf = operator[](regionI);                       \
                                                                              \
    const labelList& map =                                                    \
        mesh().typeMap(addressingTypes::FACE, regionI, regionI0);             \
                                                                              \
    forAll (vf, facei)                                                        \
    {                                                                         \
        vf0[map[facei]] = vf[facei];                                          \
    }                                                                         \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    defineRegionSurfaceFieldMapping(scalar);
    defineRegionSurfaceFieldMapping(vector);
    defineRegionSurfaceFieldMapping(tensor);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
