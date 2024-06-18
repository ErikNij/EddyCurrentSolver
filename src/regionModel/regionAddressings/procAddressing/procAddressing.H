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
    Foam::procAddressing

Description
    ...

SourceFiles
    procAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef procAddressing_H
#define procAddressing_H

#include "foamTime.H"
#include "polyMesh.H"
#include "addressingTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class procAddressing Declaration
\*---------------------------------------------------------------------------*/

class procAddressing
:
    public addressingTypes
{
private:

    // Private data

        //- Reference to time
        const Time& time_;

        //- Region name
        const word regionName_;

        //- Addressings
        mutable List<autoPtr<labelIOList> > addressingPtrs_;

        //- Face map
        mutable autoPtr<labelList> faceMapPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        procAddressing(const procAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const procAddressing&);

        // Helper functions

            //- Read addressing
            labelIOList* readAddressing(addressingType type) const;

            //- Calculate face map
            labelList* calcFaceMap() const;


public:

    //- Runtime type information
    TypeName("procAddressing");

    // Constructors

        //- Construct from mesh only
        procAddressing
        (
            const polyMesh& mesh
        );

        //- Construct from mesh and region name
        procAddressing
        (
            const polyMesh& mesh,
            const word& regionName
        );

        //- Construct from time and region name
        procAddressing
        (
            const Time& time,
            const word& regionName
        );

    // Destructor

        virtual ~procAddressing();

    // Member Functions

        //- Return const reference to time
        inline const Time& time() const
        {
            return time_;
        };

        //- Return const reference to the region name
        inline const word& regionName() const
        {
            return regionName_;
        };

        //- Return const reference to addressing of given type
        const labelList& typeAddressing(addressingType type) const;

        //- Return const reference to map of given type
        const labelList& typeMap(addressingType type) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
