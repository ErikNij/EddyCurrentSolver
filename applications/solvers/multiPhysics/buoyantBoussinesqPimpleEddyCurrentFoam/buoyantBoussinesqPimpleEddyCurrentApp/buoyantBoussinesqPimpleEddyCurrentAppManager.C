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

#include "buoyantBoussinesqPimpleEddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Settings::read() const
{
    buoyantBoussinesqPimpleEddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            buoyantBoussinesqPimpleEddyCurrentApp::Manager::debug()
        );

    buoyantBoussinesqPimpleEddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            buoyantBoussinesqPimpleEddyCurrentApp::Control::debug()
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Storage::Item_T::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::THERMAL],
        IOo
    );

    set
    (
        regionVolScalarField::LinkOrNew
        (
            IOobject
            (
                IOo.name(),
                IOo.instance(),
                IOo.db()
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimTemperature,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Storage::Item_lambda::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::THERMAL],
        IOo
    );

    set
    (
        regionVolScalarField::LinkOrNew
        (
            IOobject
            (
                IOo.name(),
                IOo.instance(),
                IOo.db()
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimThermalConductivity,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Storage::Item_rhoCp::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::THERMAL],
        IOo
    );

    set
    (
        regionVolScalarField::LinkOrNew
        (
            IOobject
            (
                IOo.name(),
                IOo.instance(),
                IOo.db()
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimDensity*dimSpecificHeatCapacity,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


// TODO: Do not write U/phi
// TODO: Substitute U with phi
void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Storage::Item_U::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::THERMAL],
        IOo
    );

    set
    (
        regionVolVectorField::LinkOrNew
        (
            IOobject
            (
                IOo.name(),
                IOo.instance(),
                IOo.db()
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVelocity,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Storage::create() const
{
    item_T().enable();
    item_lambda().enable();
    item_rhoCp().enable();
// TODO: Substitute U with phi
    item_U().enable();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Regions::create() const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Manager
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
                "Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::Manager(...) : "
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
                word(this->regionsDict().lookup("THERMAL")),
                Region::THERMAL
            );

            mesh.init(regionNameHashTable);

            this->messages().newLine();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::buoyantBoussinesqPimpleApp::Manager&
Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::buoyantBoussinesqPimpleAppManager() const
{
    if (buoyantBoussinesqPimpleAppManager_.empty())
    {
        buoyantBoussinesqPimpleAppManager_.set
        (
            new buoyantBoussinesqPimpleApp::Manager
            (
                this->args(),
                this->time(),
                this->mesh()[Region::FLUID],
                false
            )
        );
    }

    return buoyantBoussinesqPimpleAppManager_();
}

Foam::eddyCurrentApp::Manager&
Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::eddyCurrentAppManager() const
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


bool Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::setCoNum
(
    scalar& CourantNumber
) const
{
    return buoyantBoussinesqPimpleAppManager().setCoNum(CourantNumber);
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::read() const
{
    buoyantBoussinesqPimpleAppManager().read();
    eddyCurrentAppManager().read();

    // Make sure volume force is enabled
    buoyantBoussinesqPimpleAppManager().regions().region_DEFAULT().
        settings().volumeForce = true;

    // Make sure heat source is enabled
    buoyantBoussinesqPimpleAppManager().regions().region_DEFAULT().
        settings().heatSource = true;

    settings().checkRead();
    regions().checkRead();
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::init() const
{
    buoyantBoussinesqPimpleAppManager().init();
    eddyCurrentAppManager().init();

    storage().checkInit();
    regions().checkInit();
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::next() const
{
    buoyantBoussinesqPimpleAppManager().next();
    eddyCurrentAppManager().next();
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::write() const
{
    buoyantBoussinesqPimpleAppManager().write();
    eddyCurrentAppManager().write();
}


void Foam::buoyantBoussinesqPimpleEddyCurrentApp::Manager::finalize() const
{
    buoyantBoussinesqPimpleAppManager().finalize();
    eddyCurrentAppManager().finalize();
}


// ************************************************************************* //

