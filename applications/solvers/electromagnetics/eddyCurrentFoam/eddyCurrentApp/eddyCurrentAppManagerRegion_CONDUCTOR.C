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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_VRe::create() const
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

    mesh().schemesDict().setFluxRequired(name());
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_VIm::create() const
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

    mesh().schemesDict().setFluxRequired(name());
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_jRe::create() const
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

    set
    (
        new volVectorField
        (
            IOo,
            mesh(),
            dimensionedVector
            (
                word(),
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_jIm::create() const
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

    set
    (
        new volVectorField
        (
            IOo,
            mesh(),
            dimensionedVector
            (
                word(),
                dimCurrent/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_phiGradA0nRe::create() const
{
    set
    (
        new surfaceVectorField
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
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_phiGradA0nIm::create() const
{
    set
    (
        new surfaceVectorField
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
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_phiDdtARe::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimVoltage*dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_phiDdtAIm::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar
            (
                word(),
                dimVoltage*dimLength,
                0
            ),
            calculatedFvPatchScalarField::typeName
        )
    );
}


#ifdef eddyCurrentAppLink_H

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_emPrevC::create() const
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
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector
            (
                word(),
                dimLength,
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::Item_emRelDeltaA::create() const
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
                dimless,
                0.0
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_CONDUCTOR::Storage::create() const
{
    item_VRe().setState(control().meshIs3D());
    item_VIm().setState(control().meshIs3D());

    item_jRe().enable();
    item_jIm().enable();

    item_phiGradA0nRe().setState(globalSettings().biotSavart);
    item_phiGradA0nIm().setState(globalSettings().biotSavart);

    item_phiDdtARe().setState(control().meshIs3D());
    item_phiDdtAIm().setState(control().meshIs3D());

#ifdef eddyCurrentAppLink_H

    item_emPrevC().enable();
    item_emRelDeltaA().enable();

#endif
}


// ************************************************************************* //

