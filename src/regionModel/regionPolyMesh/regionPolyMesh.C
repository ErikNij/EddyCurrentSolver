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

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionPolyMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polyMesh* Foam::regionPolyMesh::newMesh(label regionI) const
{
    return new polyMesh
    (
        IOobject
        (
            regions_[regionI],
            time_.timeName(),
            time_,
            IOobject::MUST_READ
        )
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::regionPolyMesh::readRegionNames() const
{
    wordIOList regionNames
    (
        IOobject
        (
            "regions",
            this->time().constant(),
            this->time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    return regionNames;
}


void Foam::regionPolyMesh::initMeshInitialize(const wordList& regionNames) const
{
    initialized_ = false;

    regions_ = regionList(regionNames);

    meshPtrs_.clear();
    meshPtrs_.resize(regions_.size());

    addressingPtrs_.clear();
    addressingPtrs_.resize(regions_.size());
}


void Foam::regionPolyMesh::initMeshMeshes() const
{
    forAll (*this, regionI)
    {
        Info << "Create mesh for region "
            << regions()[regionI]
            << endl;

        // Create mesh
        meshPtrs_.set
        (
            regionI,
            newMesh(regionI)
        );
    }
}


void Foam::regionPolyMesh::initMeshAddressings() const
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info << "Foam::regionPolyMesh::regionPolyMesh(...) : "
                << "Create addressing for region "
                << regions()[regionI]
                << endl;
        }

        // Create addressings
        addressingPtrs_.set
        (
            regionI,
            new regionToRegionAddressing(meshPtrs_[regionI])
        );
    }
}


void Foam::regionPolyMesh::initMeshShared() const
{}


void Foam::regionPolyMesh::initMeshFinalize() const
{
    initialized_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionPolyMesh::regionPolyMesh
(
    const Time& runTime,
    bool init
)
:
    objectRegistry(
        IOobject
        (
            "regions[]",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    regions_(),
    meshPtrs_(),
    addressingPtrs_(),
    initialized_(false)
{
    if(init) this->init(readRegionNames());
}


Foam::regionPolyMesh::regionPolyMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    objectRegistry(
        IOobject
        (
            "regions[]",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    regions_(),
    meshPtrs_(),
    addressingPtrs_(),
    initialized_(false)
{
    if(init) this->init(regionNames);
}


Foam::regionPolyMesh::regionPolyMesh
(
    const Time& runTime,
    const HashTable<label>& regionNameHashTable,
    bool init
)
:
    objectRegistry(
        IOobject
        (
            "regions[]",
            runTime.constant(),
            runTime.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    time_(runTime),
    regions_(),
    meshPtrs_(),
    addressingPtrs_(),
    initialized_(false)
{
    if(init) this->init(regionNameHashTable);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionPolyMesh::~regionPolyMesh()
{
    meshPtrs_.clear();
    addressingPtrs_.clear();
}


// ************************************************************************* //

