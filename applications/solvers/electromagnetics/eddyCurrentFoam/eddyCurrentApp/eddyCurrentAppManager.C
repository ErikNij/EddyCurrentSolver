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

#include "eddyCurrentAppManager.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::eddyCurrentApp::Manager, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Settings::read() const
{
    eddyCurrentApp::Manager::debug =
        dict().lookupOrDefault
        (
            "debug",
            eddyCurrentApp::Manager::debug()
        );

    eddyCurrentApp::Control::debug =
        dict().lookupOrDefault
        (
            "debug",
            eddyCurrentApp::Control::debug()
        );

    lowFrequency = dict().lookupOrDefault("lowFrequency", false);
    biotSavart = dict().lookupOrDefault("biotSavart", false);
    lorentzForce = dict().lookupOrDefault("lorentzForce", true);
    magneticPressure = dict().lookupOrDefault("magneticPressure", true);
    jouleHeat = dict().lookupOrDefault("jouleHeat", true);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Storage::Item_sigma::create() const
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
        mesh().regions()[Region::CONDUCTOR],
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
                dimCurrent/dimVoltage/dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );

    // Make sure conductivity is not exactly zero
    scalar minSigma = VSMALL;
    if (min(get()()[Region::CONDUCTOR]).value() < minSigma)
    {
        get()()[Region::CONDUCTOR] +=
            dimensionedScalar
            (
                word(),
                dimCurrent/dimVoltage/dimLength,
                minSigma
            );
        get()()[Region::CONDUCTOR].correctBoundaryConditions();
    };
}


void Foam::eddyCurrentApp::Manager::Storage::Item_mur::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimless,
                1
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_nur::create() const
{
    set
    (
        regionVolScalarField::LinkOrNew
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                1
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_j0Re::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::DEFAULT],
        IOo
    );

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_j0Im::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", false))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::DEFAULT],
        IOo
    );

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_A0Re::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
// TODO
            fixedValueFvPatchVectorField::typeName,
//             calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_A0Im::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    if (!dict().lookupOrDefault<bool>("write", true))
    {
        IOo.writeOpt() = IOobject::NO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
// TODO
            fixedValueFvPatchVectorField::typeName,
//             calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_ARe::create() const
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
        mesh().regions()[Region::DEFAULT],
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
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_AIm::create() const
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
        mesh().regions()[Region::DEFAULT],
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
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_VReGrad::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (dict().lookupOrDefault<bool>("write", false))
    {
        IOo.writeOpt() = IOobject::AUTO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
        IOo
    );

    set
    (
        regionVolVectorField::LinkOrNew
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_VImGrad::create() const
{
    IOobject IOo
    (
        name(),
        time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (dict().lookupOrDefault<bool>("write", false))
    {
        IOo.writeOpt() = IOobject::AUTO_WRITE;
    }

    HashTable<IOobject> IOoOverride;

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
        IOo
    );

    set
    (
        regionVolVectorField::LinkOrNew
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimVoltage/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_BRe::create() const
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
        mesh().regions()[Region::DEFAULT],
        IOo
    );

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_BIm::create() const
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
        mesh().regions()[Region::DEFAULT],
        IOo
    );

    IOoOverride.set
    (
        mesh().regions()[Region::CONDUCTOR],
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
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_F::create() const
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
        mesh().regions()[Region::CONDUCTOR],
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
                dimForce/dimVolume,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_pB::create() const
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
        mesh().regions()[Region::CONDUCTOR],
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
                dimPressure,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


void Foam::eddyCurrentApp::Manager::Storage::Item_Q::create() const
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
        mesh().regions()[Region::CONDUCTOR],
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
                dimEnergy/dimTime/dimVolume,
                0
            ),
            calculatedFvPatchScalarField::typeName,
            IOoOverride
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Storage::create() const
{
    item_sigma().enable();
    item_mur().enable();
    item_nur().enable();

    item_j0Re().enable();
    item_j0Im().enable();

    item_A0Re().setState(settings().biotSavart);
    item_A0Im().setState(settings().biotSavart);

    item_ARe().enable();
    item_AIm().enable();

    item_VReGrad().setState(control().meshIs3D());
    item_VImGrad().setState(control().meshIs3D());

    item_BRe().enable();
    item_BIm().enable();

    item_F().setState(settings().lorentzForce);
    item_pB().setState(settings().magneticPressure);
    item_Q().setState(settings().jouleHeat);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Regions::create() const
{
    region_DEFAULT().enable();
    region_CONDUCTOR().enable();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Manager::Manager
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
                "Foam::eddyCurrentApp::Manager::Manager(...) : "
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

            mesh.init(regionNameHashTable);

            this->messages().newLine();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Manager::~Manager()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::next() const
{}


void Foam::eddyCurrentApp::Manager::write() const
{}


void Foam::eddyCurrentApp::Manager::finalize() const
{}


// ************************************************************************* //

