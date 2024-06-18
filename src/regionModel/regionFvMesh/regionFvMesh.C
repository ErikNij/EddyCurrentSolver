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

#include "regionFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionFvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fvMesh* Foam::regionFvMesh::newMesh(label regionI) const
{
    return new fvMesh
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


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::regionFvMesh::initMeshMeshes() const
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
            static_cast<polyMesh*>(newMesh(regionI))
        );
    }
}


void Foam::regionFvMesh::initMeshShared() const
{
    regionPolyMesh::initMeshShared();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionFvMesh::regionFvMesh
(
    const Time& runTime,
    bool init
)
:
    regionPolyMesh::regionPolyMesh
    (
        runTime,
        false
    )
{
    if(init) this->init(readRegionNames());
}


Foam::regionFvMesh::regionFvMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    regionPolyMesh::regionPolyMesh
    (
        runTime,
        false
    )
{
    if(init) this->init(regionNames);
}


Foam::regionFvMesh::regionFvMesh
(
    const Time& runTime,
    const HashTable<label>& regionNameHashTable,
    bool init
)
:
    regionPolyMesh::regionPolyMesh
    (
        runTime,
        false
    )
{
    if(init) this->init(regionNameHashTable);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionFvMesh::~regionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::regionFvMesh::operator!=(const regionFvMesh& brm) const
{
    return &brm != this;
}


bool Foam::regionFvMesh::operator==(const regionFvMesh& brm) const
{
    return &brm == this;
}


// ************************************************************************* //

