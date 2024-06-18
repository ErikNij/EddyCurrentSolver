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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interApp::Manager::Region_DEFAULT::Settings::read() const
{
    volumeForce = dict().lookupOrDefault("volumeForce", false);
    electricalConuctivity = dict().lookupOrDefault("electricalConuctivity", false);
    snGradpFromFlux = dict().lookupOrDefault("snGradpFromFlux", true);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_g::create() const
{
    set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_p::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );

    control().setpRefCell(get()());
    mesh().schemesDict().setFluxRequired(name());
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_U::create() const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_phi::create() const
{
    set
    (
        new surfaceScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(storage().U()) & mesh().Sf()
        )
    );

    get()->oldTime();
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_rho::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimDensity,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_rhoPhi::create() const
{
    set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "rho*phi",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimMass/dimTime,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_F::create() const
{
    set
    (
        new volVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimForce/dimVolume,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_alpha1::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh()
        )
    );
}

void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_alpha1Vfrac::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            storage().alpha1().weightedAverage(mesh().V())
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_transport::create() const
{
    set
    (
        new twoPhaseMixture
        (
            storage().U(),
            storage().phi(),
            storage().alpha1().name()
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_turbulence::create() const
{
    set
    (
        incompressible::turbulenceModel::New
        (
            storage().U(),
            storage().phi(),
            storage().transport()
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_interface::create() const
{
    set
    (
        new interfaceProperties
        (
            storage().alpha1(),
            storage().U(),
            storage().transport()
        )
    );
}


void Foam::interApp::Manager::Region_DEFAULT::Storage::Item_sigma::create() const
{
    set
    (
        new volScalarField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimCurrent/dimVoltage/dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_g().enable();
    item_p().enable();
    item_U().enable();
    item_phi().enable();
    item_rho().enable();
    item_rhoPhi().enable();
    item_F().setState(settings().volumeForce);

    item_alpha1().enable();
    item_alpha1Vfrac().enable();
    item_transport().enable();

    rho() = transport().rho();
    rho().correctBoundaryConditions();

    rhoPhi() = fvc::interpolate(rho()) * phi();
    rhoPhi().correctBoundaryConditions();

    item_turbulence().enable();
    item_interface().enable();

    item_sigma().setState(settings().electricalConuctivity);
}


// ************************************************************************* //

