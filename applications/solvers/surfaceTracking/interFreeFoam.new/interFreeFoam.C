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

Application
    interFreeFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "freeSurface.H"
#include "patchWave.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createDynamicFvMesh.H"

#   include "createFields.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

#       include "readPISOControls.H"

#       include "readFreeSurfaceControls.H"

        interface.moveMeshPointsForOldFreeSurfDisplacement();

        interface.updateDisplacementDirections();

//         interface.predictPoints();

        Info<< "\nMax surface Courant Number = "
            << interface.maxCourantNumber() << endl << endl;

        for (int corr=0; corr<nOuterCorr; corr++)
        {
            p.storePrevIter();

            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the fluxes relative
            phi -= fvc::meshPhi(rho, U);

#           include "CourantNo.H"

            fvVectorMatrix UEqn
            (
//                 fvm::ddt(U)
//               + fvm::div(phi, U, "div(phi,U)")
//               - fvm::laplacian(nu, U)
                fvm::ddt(rho, U)
              + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            // --- PISO loop
            for (int i=0; i<nCorr; i++)
            {
//                 volScalarField rUA = 1.0/UEqn.A();
//                 U = rUA*UEqn.H();
//                 phi = (fvc::interpolate(U) & mesh.Sf());

                volScalarField AU = UEqn.A();
                volVectorField HU = UEqn.H();

                U = HU/AU;

#               include "calcPhi.H"

// #               include "correctPhiAtInterface.H"

#               include "scalePhi.H"

                // Non-orthogonal pressure corrector loop
                for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    fvScalarMatrix pEqn
                    (
//                         fvm::laplacian(rUA, p) == fvc::div(phi)
                        fvm::laplacian(1.0/fvc::interpolate(AU), p)
                     == fvc::div(phi)
                    );

#                   include "setReference.H"

//                     pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Momentum corrector
                U -= fvc::grad(p)/AU;
//                 U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            interface.correctPoints();

#           include "freeSurfaceContinuityErrs.H"

#           include "updateSf.H"

#           include "checkResidual.H"

            AU = UEqn.A();
        }

#       include "volContinuity.H"

        Info << "Total surface tension force: "
            << interface.totalSurfaceTensionForce() << endl;

        vector totalForce =
            interface.totalViscousForce()
          + interface.totalPressureForce();

        Info << "Total force: " << totalForce << endl;

        K = interface.aMesh().faceCurvatures();

        Info << "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K)
            << ", average = " << gAverage(K) << endl << flush;

//         dist.internalField() = patchWave(mesh, patchSet, false).distance();

        runTime.write();

        if (runTime.outputTime())
        {
            interface.writeVTK();
            interface.writeVTKControlPoints();
        }

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
