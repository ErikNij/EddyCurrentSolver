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

#include "pimpleAppManager.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pimpleApp::Manager::Region_DEFAULT::Settings::read() const
{
    volumeForce = dict().lookupOrDefault("volumeForce", false);
    snGradpFromFlux = dict().lookupOrDefault("snGradpFromFlux", true);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_p::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_U::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_phi::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_transport::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_turbulence::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_F::create() const
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


void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::Item_rho::create() const
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            time().constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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
            dimensionedScalar
            (
                transportProperties.lookup("rho")
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::pimpleApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_p().enable();
    item_U().enable();
    item_phi().enable();
    item_transport().enable();
    item_turbulence().enable();
    item_F().setState(settings().volumeForce);
    item_rho().setState(settings().volumeForce);
}


// ************************************************************************* //

