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

#include "pimpleEddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::pimpleEddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pimpleEddyCurrentApp::Manager::Settings::read() const
{
    pimpleEddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            pimpleEddyCurrentApp::Manager::debug()
        );

    pimpleEddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            pimpleEddyCurrentApp::Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::pimpleEddyCurrentApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::pimpleEddyCurrentApp::Manager::Regions::create() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleEddyCurrentApp::Manager::Manager
(
    const argList& args,
    Time& time,
    regionFvMesh& mesh,
    bool master
)
:
    regionFvSolverManager
    (
        args, time, mesh, master
    )
{
    if (master)
    {
        // Init regionMesh
        if (mesh.initialized())
        {
            FatalErrorIn
            (
                "Foam::pimpleEddyCurrentApp::Manager::Manager(...) : "
            )
                << "Region mesh is already initialized. This Manager "
                << "needs to initialize the region mesh by itself. "
                << "Construct the regionMesh with 'init=false'!"
                << abort(FatalError);
        }
        else
        {
            HashTable<label> regionNameHashTable;

            regionNameHashTable.insert
            (
                polyMesh::defaultRegion,
                Region::DEFAULT
            );

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("CONDUCTOR")),
                Region::CONDUCTOR
            );

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("FLUID")),
                Region::FLUID
            );

            mesh.init(regionNameHashTable);

            this->messages().newLine();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleEddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pimpleApp::Manager&
Foam::pimpleEddyCurrentApp::Manager::pimpleAppManager() const
{
    if (pimpleAppManager_.empty())
    {
        pimpleAppManager_.set
        (
            new pimpleApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh()[Region::FLUID],
                false
            )
        );
    }

    return pimpleAppManager_();
}

Foam::eddyCurrentApp::Manager&
Foam::pimpleEddyCurrentApp::Manager::eddyCurrentAppManager() const
{
    if (eddyCurrentAppManager_.empty())
    {
        eddyCurrentAppManager_.set
        (
            new eddyCurrentApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh(),
                false
            )
        );
    }

    return eddyCurrentAppManager_();
}


bool Foam::pimpleEddyCurrentApp::Manager::setCoNum
(
    scalar& CourantNumber
) const
{
    return pimpleAppManager().setCoNum(CourantNumber);
}


void Foam::pimpleEddyCurrentApp::Manager::read() const
{
    pimpleAppManager().read();
    eddyCurrentAppManager().read();

    // Make sure volume force is enabled
    pimpleAppManager().regions().region_DEFAULT().
        settings().volumeForce = true;

    settings().checkRead();
    regions().checkRead();
}


void Foam::pimpleEddyCurrentApp::Manager::init() const
{
    pimpleAppManager().init();
    eddyCurrentAppManager().init();

    storage().checkInit();
    regions().checkInit();
}


void Foam::pimpleEddyCurrentApp::Manager::next() const
{
    pimpleAppManager().next();
    eddyCurrentAppManager().next();
}


void Foam::pimpleEddyCurrentApp::Manager::write() const
{
    pimpleAppManager().write();
    eddyCurrentAppManager().write();
}


void Foam::pimpleEddyCurrentApp::Manager::finalize() const
{
    pimpleAppManager().finalize();
    eddyCurrentAppManager().finalize();
}


// ************************************************************************* //

