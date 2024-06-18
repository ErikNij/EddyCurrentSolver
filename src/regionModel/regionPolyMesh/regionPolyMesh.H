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
    regionPolyMesh

Description
    ...

SourceFiles
    regionPolyMesh.C
    regionPolyMeshMapping.C

\*---------------------------------------------------------------------------*/

#ifndef regionPolyMesh_H
#define regionPolyMesh_H

#include "objectRegistry.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"
#include "foamTime.H"
#include "wordIOList.H"

#include "regionList.H"
#include "addressingTypes.H"
#include "regionToRegionAddressing.H"

// TODO [Idea]: Use HashPtrTable instead of PtrList???

// TODO [High]: Check toRegionI and toRegionName of regionToRegionAddressing for
//              existence!

// TODO [High]: Allow default region to be absent. WE NEED A LOT OF CHECKS!!!

// TODO [High]: Check if all region names are unique

// TODO [Low]: Remove readRegionNames()

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class regionPolyMesh
:
    public objectRegistry
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionPolyMesh(const regionPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const regionPolyMesh&);

        // Helper functions for constructor

            //- Create new mesh for region with index regionI
            polyMesh* newMesh(label regionI) const;


protected:

    // Protected data

        //- Reference to time
        const Time& time_;

        //- List of regions
        mutable regionList regions_;

        //- Mesh pointers
        mutable PtrList<polyMesh> meshPtrs_;

        //- Region to region addressings
        mutable PtrList<regionToRegionAddressing> addressingPtrs_;

        // - Is mesh initialized?
        mutable bool initialized_;

    // Protected Member Functions

        // Helper functions

// TODO: Remove readRegionNames()
            //- Read region names from file
            wordList readRegionNames() const;

            //- Init initialize
            void initMeshInitialize(const wordList& regionNames) const;

            //- Init meshes
            virtual void initMeshMeshes() const;

            //- Init addressings
            void initMeshAddressings() const;

            //- Init shared
            virtual void initMeshShared() const;

            //- Init finalize
            void initMeshFinalize() const;


public:

    // Public typedefs

        typedef polyMesh Mesh;
        typedef polyBoundaryMesh BoundaryMesh;

    //- Runtime type information
    TypeName("regionPolyMesh");

    // Constructors

        //- Construct from time
        regionPolyMesh
        (
            const Time& runTime,
            bool init = true
        );

        //- Construct from time and list of region names
        regionPolyMesh
        (
            const Time& runTime,
            const wordList& regionNames,
            bool init = true
        );

        //- Construct from time and label hash table
        regionPolyMesh
        (
            const Time& runTime,
            const HashTable<label>& regionNameHashTable,
            bool init = true
        );

    // Destructor

        virtual ~regionPolyMesh();

    // Member Functions

        // Initialization

            //- Init meshes from hash table
            inline void init(const HashTable<label>& regionNameHashTable) const
            {
                wordList regionNames = wordList(regionNameHashTable.size());

                forAllConstIter(HashTable<label>, regionNameHashTable, iter)
                {
                    regionNames[iter()] = iter.key();
                }

                init(regionNames);
            };

            //- Init meshes
            inline void init(const wordList& regionNames) const
            {
                if (!initialized())
                {
                    initMeshInitialize(regionNames);
                    initMeshMeshes();
                    initMeshAddressings();
                    initMeshShared();
                    initMeshFinalize();
                }
            };

            inline bool initialized() const
            {
                return initialized_;
            };

        // Access

            //- Return reference to this mesh
            inline const regionPolyMesh& rmesh() const
            {
                return *this;
            };

            //- Return the object registry
            inline const objectRegistry& thisDb() const
            {
                return *this;
            };

            //- Return const reference to time
            inline const Time& time() const
            {
                return time_;
            };

            //- Return const reference to global region count
            inline label size() const
            {
                return regions_.size();
            };

            //- Return const reference to regions
            inline const regionList& regions() const
            {
                return regions_;
            };

        // Active mesh data

            //- Return the vector of geometric directions in mesh
            inline const Vector<label>& geometricD
            (
                label regionI
            ) const
            {
                return operator[](regionI).geometricD();
            };

            //- Return the vector of geometric directions in mesh
            inline const Vector<label>& geometricD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]).geometricD();
            };

            //- Return the number of valid geometric dimensions in the mesh
            inline label nGeometricD
            (
                label regionI
            ) const
            {
                return operator[](regionI).nGeometricD();
            };

            //- Return the number of valid geometric dimensions in the mesh
            inline label nGeometricD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]).nGeometricD();
            };

            //- Return the vector of solved-for directions in mesh
            inline const Vector<label>& solutionD
            (
                label regionI
            ) const
            {
                return operator[](regionI).solutionD();
            };

            //- Return the vector of solved-for directions in mesh
            inline const Vector<label>& solutionD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]).solutionD();
            };

            //- Return the number of valid solved-for dimensions in the mesh
            inline label nSolutionD
            (
                label regionI
            ) const
            {
                return operator[](regionI).nSolutionD();
            };

            //- Return the number of valid solved-for dimensions in the mesh
            inline label nSolutionD
            (
                const word& regionName = polyMesh::defaultRegion
            ) const
            {
                return operator[](regions()[regionName]).nSolutionD();
            };

        // Mapping

            //- Map mesh points from default region and replace them
            //  in a given point field of non-default region
            void map
            (
                pointField& givenPoints,
                label regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region and replace them
            //  in a given point field of non-default region
            void map
            (
                pointField& givenPoints,
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region
            //  to a new point field for non-default region
            tmp<pointField> map
            (
                label regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from default region
            //  to a new point field for non-default region
            tmp<pointField> map
            (
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region and replace them
            //  in a given point field of default region
            void rmap
            (
                pointField& givenPoints,
                label regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region and replace them
            //  in a given point field of default region
            void rmap
            (
                pointField& givenPoints,
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region
            //  to a new point field for default region
            tmp<pointField> rmap
            (
                label regionI,
                const word& patchName = ""
            ) const;

            //- Map mesh points from non-default region
            //  to a new point field for default region
            tmp<pointField> rmap
            (
                const word& regionName,
                const word& patchName = ""
            ) const;

            //- Return patch mapping between two neighbour-regions
            //  based on directMappedPatch
            labelListList patchMapDirectMapped
            (
                label fromRegionI,
                label toRegionI
            ) const;

            //- Return patch mapping between two neighbour-regions
            //  based on directMappedPatch
            labelListList patchMapDirectMapped
            (
                const word& fromRegionName,
                const word& toRegionName
            ) const;

            //- Return const reference addressing of given type
            //  and given region indices
            inline const labelList& typeAddressing
            (
                addressingTypes::addressingType type,
                label fromRegionI,
                label toRegionI
            ) const
            {
                return addressingPtrs_[fromRegionI].typeAddressing
                (
                    type,
                    regions()[toRegionI]
                );
            };

            //- Return const reference addressing of given type
            //  and given region names
            inline const labelList& typeAddressing
            (
                addressingTypes::addressingType type,
                const word& fromRegionName,
                const word& toRegionName = polyMesh::defaultRegion
            ) const
            {
                return typeAddressing
                (
                    type,
                    regions()[fromRegionName],
                    regions()[toRegionName]
                );
            };

            //- Return const reference map of given type
            //  and given region indices
            inline const labelList& typeMap
            (
                addressingTypes::addressingType type,
                label fromRegionI,
                label toRegionI
            ) const
            {
                return addressingPtrs_[fromRegionI].typeMap
                (
                    type,
                    regions()[toRegionI]
                );
            };

            //- Return const reference map of given type
            //  and given region names
            inline const labelList& typeMap
            (
                addressingTypes::addressingType type,
                const word& fromRegionName,
                const word& toRegionName = polyMesh::defaultRegion
            ) const
            {
                return typeMap
                (
                    type,
                    regions()[fromRegionName],
                    regions()[toRegionName]
                );
            };

    // Member Operators

        //- Mesh access operator by index
        virtual polyMesh& operator[]
        (
            label regionI
        ) const
        {
            return meshPtrs_[regionI];
        };

        //- Mesh access operator by regionName
        inline polyMesh& operator[]
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
