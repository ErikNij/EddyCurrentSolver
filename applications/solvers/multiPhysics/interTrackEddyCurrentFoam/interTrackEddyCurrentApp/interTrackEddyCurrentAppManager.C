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

#include "interTrackEddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::interTrackEddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Manager::Settings::read() const
{
    interTrackEddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            interTrackEddyCurrentApp::Manager::debug()
        );

    interTrackEddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            interTrackEddyCurrentApp::Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Manager::Storage::create() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackEddyCurrentApp::Manager::Regions::create() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interTrackEddyCurrentApp::Manager::Manager
(
    const argList& args,
    Time& time,
    regionDynamicFvMesh& mesh,
    bool master
)
:
    regionDynamicFvSolverManager
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
                "Foam::interTrackEddyCurrentApp::Manager::Manager(...) : "
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

            regionNameHashTable.insert
            (
                word(this->regionsDict().lookup("BUFFER")),
                Region::BUFFER
            );

            mesh.init(regionNameHashTable);

            this->messages().newLine();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTrackEddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::interTrackApp::Manager&
Foam::interTrackEddyCurrentApp::Manager::interTrackAppManager() const
{
    if (interTrackAppManager_.empty())
    {
        interTrackAppManager_.set
        (
            new interTrackApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh()[Region::FLUID],
                false
            )
        );
    }

    return interTrackAppManager_();
}

Foam::eddyCurrentApp::Manager&
Foam::interTrackEddyCurrentApp::Manager::eddyCurrentAppManager() const
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


bool Foam::interTrackEddyCurrentApp::Manager::setCoNum
(
    scalar& CourantNumber
) const
{
    return interTrackAppManager().setCoNum(CourantNumber);
}


void Foam::interTrackEddyCurrentApp::Manager::read() const
{
    interTrackAppManager().read();
    eddyCurrentAppManager().read();

    // Make sure volume force is enabled
    interTrackAppManager().regions().region_DEFAULT().
        settings().volumeForce = true;

    settings().checkRead();
    regions().checkRead();
}


void Foam::interTrackEddyCurrentApp::Manager::init() const
{
    interTrackAppManager().init();
    eddyCurrentAppManager().init();

    storage().checkInit();
    regions().checkInit();
}


void Foam::interTrackEddyCurrentApp::Manager::next() const
{
    interTrackAppManager().next();
    eddyCurrentAppManager().next();
}


void Foam::interTrackEddyCurrentApp::Manager::write() const
{
    interTrackAppManager().write();
    eddyCurrentAppManager().write();
}


void Foam::interTrackEddyCurrentApp::Manager::finalize() const
{
    interTrackAppManager().finalize();
    eddyCurrentAppManager().finalize();
}


// ************************************************************************* //

