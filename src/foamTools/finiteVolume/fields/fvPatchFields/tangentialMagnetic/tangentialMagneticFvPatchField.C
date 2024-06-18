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
    tangentialMagneticFvPatchField

Description
    TODO

\*---------------------------------------------------------------------------*/

#include "tangentialMagneticFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
tangentialMagneticFvPatchField<Type>::tangentialMagneticFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    patchDirectionMixedMixedFvPatchField<Type>(p, iF)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->nHat() = this->patch().nf();
    this->normalValueFraction() = 0.0;
    this->tangentialValueFraction() = 1.0;
}


template<class Type>
tangentialMagneticFvPatchField<Type>::tangentialMagneticFvPatchField
(
    const tangentialMagneticFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    patchDirectionMixedMixedFvPatchField<Type>(ptf, p, iF, mapper)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->nHat() = this->patch().nf();
    this->normalValueFraction() = 0.0;
    this->tangentialValueFraction() = 1.0;
}


template<class Type>
tangentialMagneticFvPatchField<Type>::tangentialMagneticFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    patchDirectionMixedMixedFvPatchField<Type>
    (
        p,
        iF,
        dict,
        Field<Type>(p.size(), pTraits<Type>::zero),
        Field<Type>(p.size(), pTraits<Type>::zero),
        scalarField(p.size(), 0.0),
        scalarField(p.size(), 1.0)
    )
{}


template<class Type>
tangentialMagneticFvPatchField<Type>::tangentialMagneticFvPatchField
(
    const tangentialMagneticFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    patchDirectionMixedMixedFvPatchField<Type>(ptf, iF)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->nHat() = this->patch().nf();
    this->normalValueFraction() = 0.0;
    this->tangentialValueFraction() = 1.0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Write
template<class Type>
void tangentialMagneticFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
