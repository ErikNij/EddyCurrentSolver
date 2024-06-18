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

#include "patchEnhancedDiffusivity.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchEnhancedDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        patchEnhancedDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchEnhancedDiffusivity::patchEnhancedDiffusivity
(
    const fvMotionSolver& mSolver,
    Istream& mdData
)
:
    motionDiffusivity(mSolver),
    patchNames_(mdData),
    alpha_(readScalar(mdData)),
    basicDiffusivityPtr_(motionDiffusivity::New(mSolver, mdData))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchEnhancedDiffusivity::~patchEnhancedDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::patchEnhancedDiffusivity::operator()() const
{
    tmp<Foam::surfaceScalarField> tD = basicDiffusivityPtr_->operator()();

    surfaceScalarField& D = tD();

    const fvMesh& mesh = D.mesh();

    const polyBoundaryMesh& bdryMesh = mesh.boundaryMesh();

    forAll (patchNames_, patchNameI)
    {
        label patchID = bdryMesh.findPatchID(patchNames_[patchNameI]);

        D.boundaryField()[patchID] *= alpha_;

        const cellList& cells = mesh.cells();

        const unallocLabelList& patchCells = bdryMesh[patchID].faceCells();

        forAll(patchCells, faceI)
        {
            label curCell = patchCells[faceI];

            const labelList& curCellFaces = cells[curCell];

            forAll(curCellFaces, fI)
            {
                label curFace = curCellFaces[fI];

//                 label patchID = mesh.boundaryMesh().whichPatch(curFace);

                if(mesh.isInternalFace(curFace))
                {
                    D.internalField()[curFace] *= alpha_;
                }
            }
        }
    }

    return tD;
}


void Foam::patchEnhancedDiffusivity::correct()
{
    basicDiffusivityPtr_->correct();
}


// ************************************************************************* //
