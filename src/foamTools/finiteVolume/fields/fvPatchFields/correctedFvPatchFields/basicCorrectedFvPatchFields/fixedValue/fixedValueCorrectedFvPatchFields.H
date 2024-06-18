/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    fixedValueCorrectedFvPatchFields

Description

SourceFiles
    fixedValueCorrectedFvPatchFields.C


\*---------------------------------------------------------------------------*/

#ifndef fixedValueCorrectedFvPatchFields_H
#define fixedValueCorrectedFvPatchFields_H

#include "fixedValueCorrectedFvPatchField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef fixedValueCorrectedFvPatchField<scalar> fixedValueCorrectedFvPatchScalarField;
typedef fixedValueCorrectedFvPatchField<vector> fixedValueCorrectedFvPatchVectorField;
typedef fixedValueCorrectedFvPatchField<tensor> fixedValueCorrectedFvPatchTensorField;
typedef fixedValueCorrectedFvPatchField<sphericalTensor> fixedValueCorrectedFvPatchSphericalTensorField;
typedef fixedValueCorrectedFvPatchField<diagTensor> fixedValueCorrectedFvPatchDiagTensorField;
typedef fixedValueCorrectedFvPatchField<symmTensor> fixedValueCorrectedFvPatchSymmTensorField;
typedef fixedValueCorrectedFvPatchField<symmTensor4thOrder> fixedValueCorrectedFvPatchSymmTensor4thOrderField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
