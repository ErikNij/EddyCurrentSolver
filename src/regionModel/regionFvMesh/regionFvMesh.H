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

Class
    regionFvMesh

Description
    ...

SourceFiles
    regionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef regionFvMesh_H
#define regionFvMesh_H

#include "fvMesh.H"

#include "regionPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionFvMesh
:
    public regionPolyMesh
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionFvMesh(const regionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionFvMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            fvMesh* newMesh(label regionI) const;


protected:

    // Protected Member Functions

        // Helper functions

            //- Init meshes
            virtual void initMeshMeshes() const;

            //- Init shared
            virtual void initMeshShared() const;


public:

    // Public typedefs

        typedef fvMesh Mesh;
        typedef fvBoundaryMesh BoundaryMesh;

    //- Runtime type information
    TypeName("regionFvMesh");

    // Constructors

        //- Construct from time
        regionFvMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionFvMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

        //- Construct from time and label hash table
        regionFvMesh
        (
            const Time& runTime,
            const HashTable<label>& regionNameHashTable,
            bool init = true
        );

    // Destructor

        virtual ~regionFvMesh();

    // Member Functions

        // Access

            //- Return reference to this mesh
            inline const regionFvMesh& rmesh() const
            {
                return *this;
            };

            //- Return reference to solution dict
            inline fvSolution& solutionDict
            (
                label regionI
            ) const
            {
                return operator[](regionI).solutionDict();
            };

            //- Return reference to solution dict
            inline fvSolution& solutionDict
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]).solutionDict();
            };

    // Member Operators

        //- Mesh access operator by index
        virtual fvMesh& operator[]
        (
            label regionI
        ) const
        {
            return static_cast<fvMesh&>(meshPtrs_[regionI]);
        };

        //- Mesh access operator by regionName
        inline fvMesh& operator[]
        (
            const word& regionName
        ) const
        {
            return operator[](regions()[regionName]);
        };

        bool operator!=(const regionFvMesh&) const;
        bool operator==(const regionFvMesh&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
