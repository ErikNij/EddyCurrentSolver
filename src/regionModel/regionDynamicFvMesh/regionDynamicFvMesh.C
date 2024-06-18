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

#include "tetPointFields.H"
#include "pointFields.H"

#include "regionDynamicFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionDynamicFvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dynamicFvMesh* Foam::regionDynamicFvMesh::newMesh(label regionI) const
{
    autoPtr<dynamicFvMesh>* meshAutoPtr = new autoPtr<dynamicFvMesh>
        (
            dynamicFvMesh::New
            (
                IOobject
                (
                    regions_[regionI],
                    time_.timeName(),
                    time_,
                    IOobject::MUST_READ
                )
            )
        );

    return meshAutoPtr->ptr();
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::regionDynamicFvMesh::initMeshMeshes() const
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info << "Foam::regionDynamicFvMesh::regionDynamicFvMesh(...) : "
                << "Create mesh for region "
                << regions()[regionI]
                << endl;
        }

        // Create mesh
        meshPtrs_.set
        (
            regionI,
            static_cast<polyMesh*>(newMesh(regionI))
        );
    }
}


void Foam::regionDynamicFvMesh::initMeshShared() const
{
    regionFvMesh::initMeshShared();

    forAll (*this, regionI)
    {
        isFeMotionSolver_.resize(regions_.size(), NULL);
        isFvMotionSolver_.resize(regions_.size(), NULL);

        // Remember if type of motion solver is fe
        isFeMotionSolver_[regionI] =
            meshPtrs_[regionI].objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        // Remember if type of motion solver is fv
        isFvMotionSolver_[regionI] =
            meshPtrs_[regionI].objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionDynamicFvMesh::regionDynamicFvMesh
(
    const Time& runTime,
    bool init
)
:
    regionFvMesh::regionFvMesh
    (
        runTime,
        false
    ),
    isFeMotionSolver_(),
    isFvMotionSolver_()
{
    if(init) this->init(readRegionNames());
}


Foam::regionDynamicFvMesh::regionDynamicFvMesh
(
    const Time& runTime,
    const wordList& regionNames,
    bool init
)
:
    regionFvMesh::regionFvMesh
    (
        runTime,
        false
    ),
    isFeMotionSolver_(),
    isFvMotionSolver_()
{
    if(init) this->init(regionNames);
}


Foam::regionDynamicFvMesh::regionDynamicFvMesh
(
    const Time& runTime,
    const HashTable<label>& regionNameHashTable,
    bool init
)
:
    regionFvMesh::regionFvMesh
    (
        runTime,
        false
    ),
    isFeMotionSolver_(),
    isFvMotionSolver_()
{
    if(init) this->init(regionNameHashTable);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionDynamicFvMesh::~regionDynamicFvMesh()
{}


// ************************************************************************* //

