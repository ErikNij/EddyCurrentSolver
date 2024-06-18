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
    Foam::regionToRegionAddressing

Description
    ...

SourceFiles
    regionToRegionAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef regionToRegionAddressing_H
#define regionToRegionAddressing_H

#include "foamTime.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "HashPtrTable.H"
#include "Map.H"

#include "regionList.H"
#include "addressingTypes.H"
#include "procAddressing.H"
#include "regionAddressing.H"

// TODO: Improve implementation of self-maps (maybe put it in adressingTypes.H)

// TODO: Improve the case of (mesh().name() == regionName)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class regionToRegionAddressing Declaration
\*---------------------------------------------------------------------------*/

class regionToRegionAddressing
:
    public addressingTypes
{
private:

    // Private data

        //- Reference to time
        const Time& time_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Proc addressings
        mutable HashPtrTable<procAddressing> procAddressingPtrs_;

        //- Region addressings
        mutable HashPtrTable<regionAddressing> regionAddressingPtrs_;

        //- Addressing dictionaries and hash tables
        mutable List<HashPtrTable<labelList> > addressingPtrs_;
        mutable List<autoPtr<IOdictionary> > addressingDictPtrs_;

        //- Face map pointer hash table
        mutable HashPtrTable<labelList> faceMapPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionToRegionAddressing(const regionToRegionAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const regionToRegionAddressing&);

        // Helper functions

            //- Access proc addressing
            const procAddressing& getProcAddressing
            (
                const word& regionName
            ) const;

            //- Access region addressing
            const regionAddressing& getRegionAddressing
            (
                const word& regionName
            ) const;

            //- Calc addressing
            labelList* calcAddressing
            (
                addressingType type,
                const word& regionName
            ) const;

            //- Calc face map
            labelList* calcFaceMap
            (
                const word& regionName
            ) const;


protected:

    // Protected Member Functions


public:

    //- Runtime type information
    TypeName("regionToRegionAddressing");

    // Constructors

        //- Construct from mesh
        regionToRegionAddressing(const polyMesh& mesh);

    // Destructor

        virtual ~regionToRegionAddressing();

    // Member Functions

        //- Return const reference to time
        inline const Time& time() const
        {
            return time_;
        };

        //- Return const reference to mesh
        inline const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return const reference to addressing of given type and region name
        const labelList& typeAddressing
        (
            addressingType type,
            const word& regionName
        ) const;

        //- Return const reference to map of given type and region name
        const labelList& typeMap
        (
            addressingType type,
            const word& regionName
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
