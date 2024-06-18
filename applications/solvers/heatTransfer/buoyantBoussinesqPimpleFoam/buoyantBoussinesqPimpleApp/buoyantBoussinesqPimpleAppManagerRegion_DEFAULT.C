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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Settings::read() const
{
    volumeForce = dict().lookupOrDefault("volumeForce", false);
    heatSource = dict().lookupOrDefault("heatSource", false);
    snGradpFromFlux = dict().lookupOrDefault("snGradpFromFlux", true);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_g::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_T::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_p::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_U::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_phi::create() const
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
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_transport::create() const
{
    set
    (
        new singlePhaseTransportModel
        (
            storage().U(),
            storage().phi()
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_beta::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("beta"))
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_TRef::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("TRef"))
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_Pr::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("Pr"))
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_Prt::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("Prt"))
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_turbulence::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_alphat::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_rhok::create() const
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
            1.0 - storage().beta()*(storage().T() - storage().TRef()),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_F::create() const
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


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_rhoRef::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("rhoRef"))
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_Q::create() const
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
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimEnergy/dimTime/dimVolume,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::Item_CpRef::create() const
{
    set
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                name(),
                time().constant(),
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(storage().transport().lookup("CpRef"))
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::buoyantBoussinesqPimpleApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_g().enable();
    item_T().enable();
    item_p().enable();
    item_U().enable();
    item_phi().enable();
    item_transport().enable();
    item_beta().enable();
    item_TRef().enable();
    item_Pr().enable();
    item_Prt().enable();
    item_turbulence().enable();
    item_alphat().enable();
    item_rhok().enable();
    item_F().setState(settings().volumeForce);
    item_rhoRef().setState(settings().volumeForce || settings().heatSource);
    item_Q().setState(settings().heatSource);
    item_CpRef().setState(settings().heatSource);
}


// ************************************************************************* //

