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
    faSubMesh

Description
    Implementation of a triangulated sub-mesh of a given base finite-area
    mesh

SourceFiles
    faSubMesh.C
    faSubMeshDemandDrivenData.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef faSubMesh_H
#define faSubMesh_H

#include "fvCFD.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "faCFD.H"
#include "faMesh.H"
#include "demandDrivenData.H"
#include "addressingTypes.H"
#include "regionToRegionAddressing.H"

// TODO [High]: Points normal correction if contact angle is specified?

// TODO [High]: Curvature weighting factors?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class faSubMesh Declaration
\*---------------------------------------------------------------------------*/

class faSubMesh
{
    // Private data

        //- Base faMesh
        const polyMesh& basePolyMesh_;

        //- Base faMesh
        const faMesh& baseAreaMesh_;

        //- Sub polyMesh
        mutable polyMesh subPolyMesh_;

        //- Sub polyMesh addressing
        const regionToRegionAddressing subToBaseAddressing_;

        //- Sub faMesh
        mutable faMesh subAreaMesh_;

        //- Split points
        const pointField& subSplitPoints_;


    // Demand-driven data

        //- Sub to base face map for finite area meshes
        mutable labelList* faceSubToBaseAreaMapPtr_;

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;


    // Static Private Data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faSubMesh(const faSubMesh&);

        //- Disallow default bitwise assignment
        void operator=(const faSubMesh&);

        // Helpers

            //- Make/Collect points for sub polyMesh
            tmp<pointField> calcNewPoints() const;

        // Demand-driven data

            //- Make sub to base face map for finite area meshes
            void makeFaceSubToBaseAreaMap() const;

            //- Make face curvatures
            void makeFaceCurvatures() const;

            //- Clear geomatry related demand-driven data
            void clearGeom() const;

            //- Clear all demand-driven data
            void clearOut() const;


public:

    // Declare name of the class and it's debug switch
    TypeName("faSubMesh");


    // Constructors

        faSubMesh
        (
            const faMesh& baseAreaMesh,
            const pointField& subSplitPoints
        );


    // Destructor

        virtual ~faSubMesh();


    // Member Functions

        // Access

            //- Return access to base polyMesh
            inline const polyMesh& basePolyMesh() const
            {
                return basePolyMesh_;
            }

            //- Return access to base faMesh
            inline const faMesh& baseAreaMesh() const
            {
                return baseAreaMesh_;
            }

            //- Return access to sub polyMesh
            inline polyMesh& subPolyMesh() const
            {
                return subPolyMesh_;
            }

            //- Return access to sub faMesh
            inline faMesh& subAreaMesh() const
            {
                return subAreaMesh_;
            }

            //- Return points to split faces of base mesh
            inline const pointField& subSplitPoints() const
            {
                return subSplitPoints_;
            }

        // Mapping

            //- Sub-to-base point map
            inline const labelList& pointSubToBaseMap() const
            {
                return  subToBaseAddressing_.typeMap
                (
                    addressingTypes::POINT,
                    basePolyMesh_.name()
                );
            }

            //- Sub-to-base face map
            inline const labelList& faceSubToBaseMap() const
            {
                return  subToBaseAddressing_.typeMap
                (
                    addressingTypes::FACE,
                    basePolyMesh_.name()
                );
            }

        // Mesh motion and mophing

            //- Is sub polyMesh moving
            inline bool moving() const
            {
                return subPolyMesh().moving();
            }

            //- Move points and update faMesh
            inline bool movePoints() const
            {
                clearGeom();

                subPolyMesh().movePoints
                (
                    calcNewPoints()
                );

                return true;
            }

        // Demand-driven data

            //- Return sub to base face map for finite area meshes
            const labelList& faceSubToBaseAreaMap() const;

            //- Return face curvatures calculated on sub faMesh
            //  interpolated to base faMesh
            const areaScalarField& faceCurvatures() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
