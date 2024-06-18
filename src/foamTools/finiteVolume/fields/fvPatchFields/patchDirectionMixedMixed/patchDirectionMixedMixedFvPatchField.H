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
    patchDirectionMixedMixedFvPatchField

Description
    TODO

SourceFiles
    patchDirectionMixedMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef patchDirectionMixedMixedFvPatchField_H
#define patchDirectionMixedMixedFvPatchField_H

#include "directionMixedMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class patchDirectionMixedMixedFvPatchField
:
    public directionMixedMixedFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("patchDirectionMixedMixed");


    // Constructors

        //- Construct from patch and internal field
        patchDirectionMixedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        patchDirectionMixedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchDirectionMixedMixedFvPatchField onto
        //  a new patch
        patchDirectionMixedMixedFvPatchField
        (
            const patchDirectionMixedMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new patchDirectionMixedMixedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        patchDirectionMixedMixedFvPatchField
        (
            const patchDirectionMixedMixedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new patchDirectionMixedMixedFvPatchField<Type>(*this, iF)
            );
        }

        //- Construct from patch, internal field, dictionary and data
        patchDirectionMixedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const Field<Type>& refValue,
            const Field<Type>& refGrad,
            const scalarField& normalValueFraction,
            const scalarField& tangentialValueFraction
        );


    // Member functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "patchDirectionMixedMixedFvPatchField.C"
#else
#   ifdef xlC
#   pragma implementation("patchDirectionMixedMixedFvPatchField.C")
#   endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
