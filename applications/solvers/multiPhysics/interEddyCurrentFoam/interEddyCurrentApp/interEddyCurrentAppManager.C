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

#include "interEddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interEddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interEddyCurrentApp::Manager::Settings::read() const
{
    interEddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            interEddyCurrentApp::Manager::debug()
        );

    interEddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            interEddyCurrentApp::Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interEddyCurrentApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interEddyCurrentApp::Manager::Regions::create() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interEddyCurrentApp::Manager::Manager
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
                "Foam::interEddyCurrentApp::Manager::Manager(...) : "
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

Foam::interEddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::interApp::Manager&
Foam::interEddyCurrentApp::Manager::interAppManager() const
{
    if (interAppManager_.empty())
    {
        interAppManager_.set
        (
            new interApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh()[Region::FLUID],
                false
            )
        );
    }

    return interAppManager_();
}

Foam::eddyCurrentApp::Manager&
Foam::interEddyCurrentApp::Manager::eddyCurrentAppManager() const
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


bool Foam::interEddyCurrentApp::Manager::setCoNum
(
    scalar& CourantNumber
) const
{
    return interAppManager().setCoNum(CourantNumber);
}


void Foam::interEddyCurrentApp::Manager::read() const
{
    interAppManager().read();
    eddyCurrentAppManager().read();

    // Make sure volume force is enabled
    interAppManager().regions().region_DEFAULT().
        settings().volumeForce = true;

    // Make sure sigma is enabled
    interAppManager().regions().region_DEFAULT().
        settings().electricalConuctivity = true;

    settings().checkRead();
    regions().checkRead();
}


void Foam::interEddyCurrentApp::Manager::init() const
{
    interAppManager().init();
    eddyCurrentAppManager().init();

    storage().checkInit();
    regions().checkInit();
}


void Foam::interEddyCurrentApp::Manager::next() const
{
    interAppManager().next();
    eddyCurrentAppManager().next();
}


void Foam::interEddyCurrentApp::Manager::write() const
{
    interAppManager().write();
    eddyCurrentAppManager().write();
}


void Foam::interEddyCurrentApp::Manager::finalize() const
{
    interAppManager().finalize();
    eddyCurrentAppManager().finalize();
}


// ************************************************************************* //

