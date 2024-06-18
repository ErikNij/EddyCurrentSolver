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

#include "interAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interApp::Manager::Settings::read() const
{
    interApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            interApp::Manager::debug()
        );

    Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interApp::Manager::Regions::create() const
{
    region_DEFAULT().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::interApp::Manager::setCoNum(scalar& CourantNumber) const
{
    CourantNumber = 0.0;

    const Time& runTime = time();
    const fvMesh& mesh = this->mesh();

    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

    tmp<surfaceScalarField> tphi(storage.phi());
    surfaceScalarField& phi = tphi();

    // Convective Courant Number
    {
#       include "CourantNo.H"

        CourantNumber = max(CourantNumber, CoNum);
    }

    tphi.clear();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interApp::Manager::Manager
(
    const argList& args,
    Time& time,
    fvMesh& mesh,
    bool master
)
:
    fvSolverManager
    (
        args, time, mesh, master
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interApp::Manager::next() const
{}


void Foam::interApp::Manager::write() const
{
    const fvMesh& mesh = this->mesh();

    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

// TODO: Remove after debug

    volVectorField gradAlpha("gradAlpha", fvc::grad(storage.alpha1()));
    gradAlpha.write();

    volScalarField magGradAlpha("magGradAlpha", mag(gradAlpha));
    magGradAlpha /= max(magGradAlpha)
        + dimensionedScalar(word(), magGradAlpha.dimensions(), VSMALL);
    magGradAlpha.write();

    volScalarField magSqrGradAlpha("magSqrGradAlpha", magSqr(gradAlpha));
    magSqrGradAlpha /= max(magSqrGradAlpha)
        + dimensionedScalar(word(), magSqrGradAlpha.dimensions(), VSMALL);
    magSqrGradAlpha.write();

    volScalarField UGradAlpha("UGradAlpha", mag(gradAlpha & storage.U()));
    UGradAlpha /= max(UGradAlpha)
        + dimensionedScalar(word(), UGradAlpha.dimensions(), VSMALL);
    UGradAlpha.write();

//     volScalarField cutFct("cutFct", 1.0-magSqrGradAlpha);
    volScalarField cutFct("cutFct", 1.0-pow(magSqrGradAlpha, 2));
//     volScalarField cutFct("cutFct", 1.0-pow(magSqrGradAlpha, 4));
    cutFct /= max(cutFct)
        + dimensionedScalar(word(), cutFct.dimensions(), VSMALL);
    cutFct.write();

    volScalarField cutFctAverage("cutFctAverage", fvc::average(fvc::interpolate(cutFct)));
    cutFctAverage.write();

    storage.rho().write();

    if (regions().region_DEFAULT().settings().electricalConuctivity)
    {
        storage.sigma().write();
    };

    storage.interface().K().write();

    volVectorField surfaceTensionForce
    (
        "FS",
        fvc::reconstruct
        (
            fvc::interpolate(storage.interface().sigmaK())
          * fvc::snGrad(storage.alpha1())
          * mesh.magSf()
         )
    );
    surfaceTensionForce.write();
}


void Foam::interApp::Manager::finalize() const
{}


// ************************************************************************* //

