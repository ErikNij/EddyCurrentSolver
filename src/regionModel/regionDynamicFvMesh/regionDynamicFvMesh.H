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
    regionDynamicFvMesh

Description
    ...

SourceFiles
    regionDynamicFvMesh.C
    regionDynamicFvMeshMapping.C

\*---------------------------------------------------------------------------*/

#ifndef regionDynamicFvMesh_H
#define regionDynamicFvMesh_H

#include "dynamicFvMesh.H"

#include "regionFvMesh.H"

// TODO [Low]: What happens if neighbouring regions do see different motion
//             solver types?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class regionDynamicFvMesh
:
    public regionFvMesh
{
private:

    // Private data

        //- Finite element motion solver type
        mutable boolList isFeMotionSolver_;

        //- Finite volume motion solver type
        mutable boolList isFvMotionSolver_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionDynamicFvMesh(const regionDynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionDynamicFvMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            dynamicFvMesh* newMesh(label regionI) const;


protected:

    // Protected Member Functions

        // Helper functions

            //- Init meshes
            virtual void initMeshMeshes() const;

            //- Init shared
            virtual void initMeshShared() const;


public:

    //- Runtime type information
    TypeName("regionDynamicFvMesh");

    // Constructors

        //- Construct from time (without init)
        regionDynamicFvMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionDynamicFvMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

        //- Construct from time and label hash table
        regionDynamicFvMesh
        (
            const Time& runTime,
            const HashTable<label>& regionNameHashTable,
            bool init = true
        );

    // Destructor

        virtual ~regionDynamicFvMesh();

    // Member Functions

        // Access

            //- Return reference to this mesh
            inline const regionDynamicFvMesh& rmesh() const
            {
                return *this;
            };

            //- Return true if motion solver type is fe
            inline bool isFeMotionSolver
            (
                label regionI
            ) const
            {
                return isFeMotionSolver_[regionI];
            };

            //- Return true if motion solver type is fe
            inline bool isFeMotionSolver
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return isFeMotionSolver_[regions()[regionName]];
            };

            //- Return true if motion solver type is fv
            inline bool isFvMotionSolver
            (
                label regionI
            ) const
            {
                return isFvMotionSolver_[regionI];
            };

            //- Return true if motion solver type is fv
            inline bool isFvMotionSolver
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return isFvMotionSolver_[regions()[regionName]];
            };

        // Mapping

            //- Map boundary field of mesh velocity between
            //  two neighbour-regions based on directMappedPatch
            void patchMapMeshVelocityDirectMapped
            (
                label fromRegionI,
                label toRegionI
            ) const;

            //- Map boundary field of mesh velocity between
            //  two neighbour-regions based on directMappedPatch
            void patchMapMeshVelocityDirectMapped
            (
                const word& fromRegionName,
                const word& toRegionName
            ) const;

            //- TODO
            void patchMapMeshPointsDirectMapped
            (
                label fromRegionI,
                label toRegionI,
                pointField& toPoints
            ) const;

            //- TODO
            void patchMapMeshPointsDirectMapped
            (
                const word& fromRegionName,
                const word& toRegionName,
                pointField& toPoints
            ) const;

    // Member Operators

        //- Mesh access operator by index
        virtual dynamicFvMesh& operator[]
        (
            label regionI
        ) const
        {
            return static_cast<dynamicFvMesh&>(meshPtrs_[regionI]);
        };

        //- Mesh access operator by regionName
        inline dynamicFvMesh& operator[]
        (
            const word& regionName
        ) const
        {
            return operator[](regions()[regionName]);
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
