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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Region_DEFAULT::Settings::read() const
{
    volumeForce = dict().lookupOrDefault("volumeForce", false);
    snGradpFromFlux = dict().lookupOrDefault("snGradpFromFlux", true);
    cTransport = dict().lookupOrDefault("cTransport", false);
    relToUinf = dict().lookupOrDefault("relToUinf", false);
    creepingFlow = dict().lookupOrDefault("creepingFlow", false);
    heleShawPoissonDrag = dict().lookupOrDefault("heleShawPoissonDrag", false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_g::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_p::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_U::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_phi::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_rho::create() const
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
                IOobject::AUTO_WRITE
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_F::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_fluidIndicator::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_transport::create() const
{
    set
    (
        new twoPhaseMixture
        (
            storage().U(),
            storage().phi(),
            storage().fluidIndicator().name()
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_turbulence::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_interface::create() const
{
// TODO: Add more constructors and simplify
//       Make it use real rho and mu fields!
    set
    (
        new trackedSurface
        (
            storage().mesh(),
            storage().rho(),
            storage().U(),
            storage().p(),
            storage().phi(),
            NULL,
            storage().item_c().getPtr(),
            storage().item_g().getPtr(),
            storage().item_transport().getPtr(),
            storage().item_turbulence().getPtr(),
            NULL
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_Dc::create() const
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
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_c::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_Uinf::create() const
{
    set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_UinfOld::create() const
{
    set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedVector
            (
                word(),
                dimVelocity,
                vector::zero
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_fInf::create() const
{
    set
    (
        new uniformDimensionedVectorField
        (
            IOobject
            (
                name(),
                time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedVector
            (
                word(),
                dimForce,
                vector::zero
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_gradcInf::create() const
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


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_dUinfRelax::create() const
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
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_mInf::create() const
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
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::Item_heleShawGapWidth::create() const
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
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interTrackApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_g().enable();
    item_p().enable();
    item_U().enable();
    item_phi().enable();
    item_rho().enable();
    item_F().setState(settings().volumeForce);

    item_fluidIndicator().enable();
    item_transport().enable();

    rho() = transport().rho();
    rho().correctBoundaryConditions();

    item_turbulence().enable();

    item_Dc().setState(settings().cTransport);
    item_c().setState(settings().cTransport);

    item_Uinf().setState(settings().relToUinf);
    item_UinfOld().setState(settings().relToUinf);
    item_fInf().setState(settings().relToUinf);
    item_gradcInf().setState(settings().relToUinf);
    item_dUinfRelax().setState(settings().relToUinf);
    item_mInf().setState(settings().relToUinf);

    item_heleShawGapWidth().setState(settings().heleShawPoissonDrag);

    item_interface().enable();
}


// ************************************************************************* //

