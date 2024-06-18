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

Description
    Owner/Neighbour "interpolation" static functions

\*---------------------------------------------------------------------------*/

#ifndef ownerNeighbour_H
#define ownerNeighbour_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
ownerInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        surfaceScalarField
        (
            IOobject
            (
                "ownerWeightingFactors" + vf.name(),
                vf.mesh().time().timeName(),
                vf.mesh()
            ),
            vf.mesh(),
            dimensionedScalar(word(), dimless, 1.0)
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
ownerInterpolate
(const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
        ownerInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
neighbourInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        surfaceScalarField
        (
            IOobject
            (
                "neighbourWeightingFactors" + vf.name(),
                vf.mesh().time().timeName(),
                vf.mesh()
            ),
            vf.mesh(),
            dimensionedScalar(word(), dimless, 0.0)
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
neighbourInterpolate
(const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
        neighbourInterpolate(tvf());
    tvf.clear();
    return tinterp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
