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


#include "freeSurface.H"

#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void freeSurface::moveFixedPatches(const vectorField& displacement)
{
    // Take only displacement at fixed patches
    vectorField delta(aMesh().nPoints(), vector::zero);

    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID =
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("freeSurface::moveFixedPatches(...)")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            delta[patchPoints[pointI]] = displacement[patchPoints[pointI]];
        }
    }


    //-- Set mesh motion boundary conditions

    // Check mesh motion solver type

    bool feMotionSolver =
        mesh().objectRegistry::foundObject<tetPointVectorField>
        (
            "motionU"
        );

    bool fvMotionSolver =
        mesh().objectRegistry::foundObject<pointVectorField>
        (
            "pointMotionU"
        );

    if (feMotionSolver)
    {
        tetPointVectorField& motionU =
            const_cast<tetPointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<tetPointVectorField>
                (
                    "motionU"
                )
            );

        fixedValueTetPolyPatchVectorField& motionUaPatch =
            refCast<fixedValueTetPolyPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        tetPolyPatchInterpolation tppiAPatch
        (
            refCast<const faceTetPolyPatch>
            (
                motionUaPatch.patch()
            )
        );

        motionUaPatch ==
            tppiAPatch.pointToPointInterpolate
            (
                delta/DB().deltaT().value()
            );

        if(twoFluids_)
        {
            pointField deltaB =
                interpolatorAB().pointInterpolate
                (
                    delta
                );

            fixedValueTetPolyPatchVectorField& motionUbPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            tetPolyPatchInterpolation tppiBPatch
            (
                refCast<const faceTetPolyPatch>(motionUbPatch.patch())
            );

            motionUbPatch ==
                tppiBPatch.pointToPointInterpolate
                (
                    deltaB/DB().deltaT().value()
                );
        }
    }
    else if (fvMotionSolver)
    {
        pointVectorField& motionU =
            const_cast<pointVectorField&>
            (
                mesh().objectRegistry::
                lookupObject<pointVectorField>
                (
                    "pointMotionU"
                )
            );

        fixedValuePointPatchVectorField& motionUaPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                motionU.boundaryField()[aPatchID()]
            );

        motionUaPatch ==
            delta/DB().deltaT().value();

        if(twoFluids_)
        {
            pointField deltaB =
                interpolatorAB().pointInterpolate
                (
                    displacement
                );

            fixedValuePointPatchVectorField& motionUbPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[bPatchID()]
                );

            motionUbPatch ==
                deltaB/DB().deltaT().value();
        }
    }

    mesh().update();

    aMesh().movePoints();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
