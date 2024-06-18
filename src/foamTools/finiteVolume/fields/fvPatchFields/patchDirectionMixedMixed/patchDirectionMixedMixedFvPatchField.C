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

\*---------------------------------------------------------------------------*/

#include "patchDirectionMixedMixedFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
patchDirectionMixedMixedFvPatchField<Type>::patchDirectionMixedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    directionMixedMixedFvPatchField<Type>(p, iF)
{
    this->nHat() = this->patch().nf();
}


template<class Type>
patchDirectionMixedMixedFvPatchField<Type>::patchDirectionMixedMixedFvPatchField
(
    const patchDirectionMixedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedMixedFvPatchField<Type>(ptf, p, iF, mapper)
{
    this->nHat() = this->patch().nf();
}


template<class Type>
patchDirectionMixedMixedFvPatchField<Type>::patchDirectionMixedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedMixedFvPatchField<Type>
    (
        p,
        iF,
        dict,
        Field<Type>("refValue", dict, p.size()),
        Field<Type>("refGradient", dict, p.size()),
        vectorField(p.nf()),
        scalarField("normalValueFraction", dict, p.size()),
        scalarField("tangentialValueFraction", dict, p.size())
    )
{}


template<class Type>
patchDirectionMixedMixedFvPatchField<Type>::patchDirectionMixedMixedFvPatchField
(
    const patchDirectionMixedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    directionMixedMixedFvPatchField<Type>(ptf, iF)
{
    this->nHat() = this->patch().nf();
}


template<class Type>
patchDirectionMixedMixedFvPatchField<Type>::patchDirectionMixedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const Field<Type>& refValue,
    const Field<Type>& refGrad,
    const scalarField& normalValueFraction,
    const scalarField& tangentialValueFraction
)
:
    directionMixedMixedFvPatchField<Type>
    (
        p,
        iF,
        dict,
        Field<Type>(refValue),
        Field<Type>(refGrad),
        vectorField(p.nf()),
        scalarField(normalValueFraction),
        scalarField(tangentialValueFraction)
    )
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Write
template<class Type>
void patchDirectionMixedMixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->refValue().writeEntry("refValue", os);
    this->refGrad().writeEntry("refGradient", os);
    this->normalValueFraction().writeEntry("normalValueFraction", os);
    this->tangentialValueFraction().writeEntry("tangentialValueFraction", os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
