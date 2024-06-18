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

#include "calculatedFvPatchFields.H"
#include "extrapolatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void
extrapolate
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    // Create list with boundary types
    wordList vfPatchTypes(mesh.boundary().size());

    forAll(vf.boundaryField(), patchI)
    {
        if
        (
            vf.boundaryField()[patchI].type()
         == calculatedFvPatchField<Type>::typeName
        )
        {
            vfPatchTypes[patchI] =
                extrapolatedFvPatchField<Type>::typeName;
        }
        else
        {
            vfPatchTypes[patchI] =
                vf.boundaryField()[patchI].type();
        }
    }

    // Create copy of vf with extrapolated patches
    tmp< GeometricField<Type, fvPatchField, volMesh> > tvfEx
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                vf.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            vf,
            vfPatchTypes
        )
    );

    // Extrapolation
    tvfEx().correctBoundaryConditions();

    // Copy extrapolated field and delete tmp
    vf == tvfEx;
}

template<class Type>
void
extrapolate
(
    tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    extrapolate(vf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
