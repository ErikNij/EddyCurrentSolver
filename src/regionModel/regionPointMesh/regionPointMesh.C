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

#include "regionPointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionPointMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointMesh* Foam::regionPointMesh::newMesh(label regionI) const
{
    return new pointMesh(rpMesh_[regionI]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionPointMesh::regionPointMesh
(
    const regionPolyMesh& rpMesh
)
:
    rpMesh_(rpMesh),
    meshPtrs_()
{
    forAll (*this, regionI)
    {
        if (debug)
        {
            Info << "Foam::regionPointMesh::regionPointMesh(...) : "
                << "Create mesh for region "
                << regions()[regionI]
                << endl;
        }

        // Create mesh
        meshPtrs_.set
        (
            regionI,
            newMesh(regionI)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionPointMesh::~regionPointMesh()
{
    meshPtrs_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::regionPointMesh::operator!=(const regionPointMesh& brm) const
{
    return &brm != this;
}


bool Foam::regionPointMesh::operator==(const regionPointMesh& brm) const
{
    return &brm == this;
}


// ************************************************************************* //

