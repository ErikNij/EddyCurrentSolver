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

#include "relRefDeltaDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relRefDeltaDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        relRefDeltaDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relRefDeltaDiffusivity::relRefDeltaDiffusivity
(
    const fvMotionSolver& mSolver,
    Istream& mdData
)
:
    uniformDiffusivity(mSolver, mdData),
    refDelta_
    (
        IOobject
        (
            "refDelta",
            mSolver.mesh().time().constant(),
            mSolver.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        1.0/mSolver.mesh().deltaCoeffs()
    )
{
    if (!refDelta_.headerOk())
    {
        refDelta_.write();
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relRefDeltaDiffusivity::~relRefDeltaDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::relRefDeltaDiffusivity::correct()
{
    const fvMesh& mesh = mSolver().mesh();

    surfaceScalarField delta
    (
        IOobject
        (
            "delta",
            mesh.time().timeName(),
            mesh
        ),
        1.0/mesh.deltaCoeffs()
    );

    surfaceScalarField relDiffDelta = mag((delta-refDelta_)/refDelta_);

    relDiffDelta += dimensionedScalar(word(), dimless, SMALL);

    faceDiffusivity_ = relDiffDelta;
}


// ************************************************************************* //
