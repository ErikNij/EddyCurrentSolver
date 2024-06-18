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

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Settings::read() const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_f0::create() const
{
    if (globalSettings().freqOpt)
    {
        set
        (
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    name(),
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                dimensionedScalar
                (
                    word(),
                    dimless/dimTime,
                    0.0
                )
            )
        );
    }
    else
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
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_omega0::create() const
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
            mathematicalConstant::twoPi*storage().f0()
        )
    );
}


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_phiGradAnRe::create() const
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


void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_phiGradAnIm::create() const
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


// TEST
// void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_GRe::create() const
// {
//     set
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 name(),
//                 time().timeName(),
//                 mesh(),
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh()
//         )
//     );
// }
//
//
// void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::Item_GIm::create() const
// {
//     set
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 name(),
//                 time().timeName(),
//                 mesh(),
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh()
//         )
//     );
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::eddyCurrentApp::Manager::Region_DEFAULT::Storage::create() const
{
    item_f0().enable();
    item_omega0().enable();

    item_phiGradAnRe().enable();
    item_phiGradAnIm().enable();

// TEST
//     item_GRe().enable();
//     item_GIm().enable();
}


// ************************************************************************* //

