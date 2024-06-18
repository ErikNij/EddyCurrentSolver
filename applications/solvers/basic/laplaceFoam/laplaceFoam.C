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
    laplacianFoam

Description
    Solves Laplacian.

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readSIMPLEControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (!args.optionFound("overwrite"))
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
    }

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix psiEqn
        (
            fvm::laplacian(gamma, psi) == cellsource
        );

        if (!deflated)
        {
            psiEqn.setReference(psiRefCell, psiRefValue);
        }

        psiEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi = psiEqn.flux();
        }
    }

    psiGrad == fvc::grad(psi);
    psiGradFlux == fvc::reconstruct(phi/gamma);
    psiGradSnGrad == fvc::reconstruct(fvc::snGrad(psi) * mesh.magSf());

    // Write
    runTime.writeNow();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
