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

#include "interTrackAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interTrackApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Settings::read() const
{
    interTrackApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            interTrackApp::Manager::debug()
        );

    Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Regions::create() const
{
    region_DEFAULT().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Something seems to be wrong here (after bulk mesh movement)!
bool Foam::interTrackApp::Manager::setCoNum(scalar& CourantNumber) const
{
    CourantNumber = 0.0;

    const Time& runTime = time();
    const fvMesh& mesh = this->mesh();

    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

    tmp<surfaceScalarField> tphi(storage.phi());
    surfaceScalarField& phi = tphi();

    // Convective Courant Number
    {
        if (mesh.moving())
        {
            const volVectorField& U = storage.U();

// TODO: Something is wrong with meshPhi
            // Make fluxes relative
            phi -= fvc::meshPhi(U);
// // TODO: Use U to calculate phi for now. Something with meshPhi is totally wrong
//             // Calculate phi from U
//             phi = fvc::interpolate(storage.U()) & mesh.Sf();

        }

#       include "CourantNo.H"

        CourantNumber = max(CourantNumber, CoNum);
    }

    // Mesh Courant Number
    if (mesh.moving())
    {
#       include "meshCourantNo.H"

        CourantNumber = max(CourantNumber, meshCoNum);
    }

    // Interface Courant Number
    {
// TODO: Use const
        trackedSurface& interface = storage.interface();

        if (!interface.fixedInterface())
        {
            scalar interfaceCoNum = interface.maxCourantNumber();

            Info << "Surface Courant Number max: "
                << interface.maxCourantNumber() << endl;

            CourantNumber = max(CourantNumber, interfaceCoNum);
        }
    }

    tphi.clear();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interTrackApp::Manager::Manager
(
    const argList& args,
    Time& time,
    dynamicFvMesh& mesh,
    bool master
)
:
    dynamicFvSolverManager
    (
        args, time, mesh, master
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTrackApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::next() const
{
    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

// TODO TEST: Sub-mesh
#   include "trackedSurfaceStatistics.H"
}


void Foam::interTrackApp::Manager::write() const
{
    Region_DEFAULT::Storage& storage = regions().region_DEFAULT().storage();

// TODO: Use const
    trackedSurface& interface = storage.interface();

    if (debug)
    {
        interface.writeVTK();
        interface.writeA();
        interface.writeVolA();
    }
}


void Foam::interTrackApp::Manager::finalize() const
{}


// ************************************************************************* //

