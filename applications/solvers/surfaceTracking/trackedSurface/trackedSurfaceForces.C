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

#include "trackedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

vectorField trackedSurface::pressureForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals();

    const scalarField& P = p().boundaryField()[aPatchID()];

//     vectorField pressureForce = S * n & (P*I);

    vectorField pressureForce = S * n * P;

    return pressureForce;
}


vectorField trackedSurface::viscousForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals();

    symmTensorField devReff =
        rho().boundaryField()[aPatchID()]
      * turbulence().devReff()().boundaryField()[aPatchID()];

    vectorField viscousForce = S * n & devReff;

//     vectorField nGradU = U().boundaryField()[aPatchID()].snGrad();
//
//     vectorField viscousForce = -nGradU - n*nGradUn();
//
//     if (!fixedInterface_)
//     {
//          viscousForce -= (fac::grad(Us())().internalField()&n);
//     }
//
//     viscousForce *= muEffFluidAval() * S;

    return viscousForce;
}

vector trackedSurface::totalPressureForce() const
{
    return gSum(pressureForce());
}


vector trackedSurface::totalViscousForce() const
{
    return gSum(viscousForce());
}


vector trackedSurface::normalViscousForce() const
{
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField force = viscousForce();

    vectorField normalViscousForce = (sqr(n)&force);

    return gSum(normalViscousForce);
}


vector trackedSurface::tangentialViscousForce() const
{
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField force = viscousForce();

    vectorField tangentialViscousForce = ((I-sqr(n))&force);

    return gSum(tangentialViscousForce);
}


vector trackedSurface::totalNormalForce() const
{
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField force = pressureForce() + viscousForce();

    vectorField normalForce = (sqr(n)&force);

    return gSum(normalForce);
}


vector trackedSurface::totalTangentialForce() const
{
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField force = pressureForce() + viscousForce();

    vectorField tangentialForce = ((I-sqr(n))&force);

    return gSum(tangentialForce);
}


vector trackedSurface::totalForce() const
{
    vectorField force = pressureForce() + viscousForce();

    return gSum(force);
}


vector trackedSurface::totalNormalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField normalSurfaceTensionForce =
        S*(sqr(n)&surfaceTensionForce());

    return gSum(normalSurfaceTensionForce);
}


vector trackedSurface::totalTangentialSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals();

    vectorField tangentialSurfaceTensionForce =
        S*((I-sqr(n))&surfaceTensionForce());

    return gSum(tangentialSurfaceTensionForce);
}


vector trackedSurface::totalSurfaceTensionForce() const
{
    return gSum(surfaceTensionForce());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
