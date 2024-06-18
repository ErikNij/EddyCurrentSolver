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

#include "buoyantBoussinesqPimpleAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::buoyantBoussinesqPimpleApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Settings::read() const
{
    buoyantBoussinesqPimpleApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            buoyantBoussinesqPimpleApp::Manager::debug()
        );

    Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Regions::create() const
{
    region_DEFAULT().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::buoyantBoussinesqPimpleApp::Manager::setCoNum(scalar& CourantNumber) const
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

Foam::buoyantBoussinesqPimpleApp::Manager::Manager
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

Foam::buoyantBoussinesqPimpleApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::next() const
{}


void Foam::buoyantBoussinesqPimpleApp::Manager::write() const
{}


void Foam::buoyantBoussinesqPimpleApp::Manager::finalize() const
{}


// ************************************************************************* //

