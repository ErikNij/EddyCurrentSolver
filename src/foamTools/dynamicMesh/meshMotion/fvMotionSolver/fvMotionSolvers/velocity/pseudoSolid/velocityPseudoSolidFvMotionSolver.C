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

#include "velocityPseudoSolidFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"
#include "leastSquaresVolPointInterpolation.H"
#include "fvc.H"
#include "slipFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityPseudoSolidFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        velocityPseudoSolidFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityPseudoSolidFvMotionSolver::velocityPseudoSolidFvMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    fvMotionSolver(mesh),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(fvMesh_)
    ),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            vector::zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(*this, lookup("diffusivity"))
    ),
    leastSquaresVolPoint_(false)
{
    if (found("leastSquaresVolPoint"))
    {
        leastSquaresVolPoint_ =
            Switch(lookup("leastSquaresVolPoint"));
    }

    const dictionary& pseudoSolidDic = subDict("pseudoSolid");

    nu_ = readScalar(pseudoSolidDic.lookup("poissonsRatio"));

    nCorrectors_ =  readInt(pseudoSolidDic.lookup("nCorrectors"));

    convergenceTolerance_ =
        readScalar(pseudoSolidDic.lookup("convergenceTolerance"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityPseudoSolidFvMotionSolver::~velocityPseudoSolidFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityPseudoSolidFvMotionSolver::curPoints() const
{
    if (leastSquaresVolPoint_)
    {
        leastSquaresVolPointInterpolation::New(fvMesh_).interpolate
        (
            cellMotionU_,
            pointMotionU_
        );
    }
    else
    {
        volPointInterpolation::New(fvMesh_).interpolate
        (
            cellMotionU_,
            pointMotionU_
        );
    }

    tmp<pointField> tcurPoints(new pointField(fvMesh_.allPoints()));
    pointField& cp = tcurPoints();
    const pointField& pointMotionUI = pointMotionU_.internalField();

    forAll(pointMotionUI, pointI)
    {
        cp[pointI] +=
            pointMotionUI[pointI]*fvMesh_.time().deltaT().value();
    }

    // tmp<pointField> tcurPoints
    // (
    //     fvMesh_.points()
    //   + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
    // );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::velocityPseudoSolidFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.allPoints());

    diffusivityPtr_->correct();

    // ZT, Problem on symmetry plane
//     pointMotionU_.boundaryField().updateCoeffs();

    int iCorr = 0;
    scalar initialResidual = 0;

    do
    {
        Info << "Correction: " << ++iCorr << endl;

        surfaceScalarField muf = diffusivityPtr_->operator()()();
        surfaceScalarField mutf = muf;
        surfaceScalarField lambdaf(word(), mutf*(2*nu_/(1 - 2*nu_)));

        scalar bs = 0.1;
        const unallocLabelList& owner = fvMesh_.owner();
        forAll(owner, faceI)
        {
            mutf[faceI] *= bs;
            lambdaf[faceI] *= bs;
        }

        volScalarField mu
        (
            "mu",
            fvc::average(muf)
          / dimensionedScalar(word(), muf.dimensions(), 1.0)
        );

        volScalarField mut
        (
            "mut",
            fvc::average(mutf)
          / dimensionedScalar(word(), muf.dimensions(), 1.0)
        );

        volScalarField lambda
        (
            "lambda",
            fvc::average(lambdaf)
          / dimensionedScalar(word(), lambdaf.dimensions(), 1.0)
        );

        volVectorField& U = cellMotionU_;
        volTensorField gradU("gradCellMotionU", fvc::grad(U));

        fvVectorMatrix motionEqn
        (
            fvm::laplacian
            (
                mu + (mut + lambda), U,
                "laplacian(diffusivity,cellMotionU)"
            )
          + fvc::div
            (
                mut*gradU.T() + lambda*(I*tr(gradU)) - (mut + lambda)*gradU,
                "div(cellMotionSigma)"
            )
        );

        // Solve the motion equation
        initialResidual = motionEqn.solve().initialResidual();

        Info << "Initial residual: " << initialResidual << endl;
    }
    while (initialResidual > convergenceTolerance_ && iCorr < nCorrectors_);
}


void Foam::velocityPseudoSolidFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    fvMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(NULL);
    diffusivityPtr_ = motionDiffusivity::New(*this, lookup("diffusivity"));
}


// ************************************************************************* //
