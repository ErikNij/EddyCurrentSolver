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

#include "relRefVolumeDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relRefVolumeDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        relRefVolumeDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relRefVolumeDiffusivity::relRefVolumeDiffusivity
(
    const fvMotionSolver& mSolver,
    Istream& mdData
)
:
    uniformDiffusivity(mSolver, mdData),
    refV_
    (
        IOobject
        (
            "refV",
            mSolver.mesh().time().constant(),
            mSolver.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mSolver.mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    )
{
    if (!refV_.headerOk())
    {
        refV_.internalField() = mSolver.mesh().V();
        refV_.correctBoundaryConditions();

        refV_.write();
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relRefVolumeDiffusivity::~relRefVolumeDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::relRefVolumeDiffusivity::correct()
{
    const fvMesh& mesh = mSolver().mesh();

    volScalarField V
    (
        IOobject
        (
            "delta",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );

    V.internalField() = mesh.V();
    V.correctBoundaryConditions();

    volScalarField relDiffV = mag((V-refV_)/refV_);

    relDiffV += dimensionedScalar(word(), dimless, SMALL);

    faceDiffusivity_ = fvc::interpolate(relDiffV, "interpolate(V)");
}


// ************************************************************************* //
