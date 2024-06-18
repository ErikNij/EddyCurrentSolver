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
    Foam::addressingTypes

Description
    ...

SourceFiles
    addressingTypes.C

\*---------------------------------------------------------------------------*/

#ifndef addressingTypes_H
#define addressingTypes_H

#include "NamedEnum.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class addressingTypes Declaration
\*---------------------------------------------------------------------------*/

class addressingTypes
{
public:

    //- Addressing type enums
    enum addressingType
    {
        POINT,
        FACE,
        CELL,
        BOUNDARY
    };

    // Static data

        static const NamedEnum<addressingType, 4> addressingNames;


protected:

    // Protected static inline member funcions

        static inline label size
        (
            addressingType type,
            const polyMesh& mesh
        )
        {
            switch ( type )
            {
                case POINT: return mesh.points().size();
                case FACE: return mesh.faces().size();
                case CELL: return mesh.cells().size();
                case BOUNDARY: return mesh.boundaryMesh().size();
                default: return -1;
            }
        }

        static inline label invalidAddressingLabel
        (
            addressingType type
        )
        {
            return (type == FACE) ? 0 : -1;
        }

        static inline label invalidAddressing
        (
            addressingType type,
            label typeI
        )
        {
            return typeI == invalidAddressingLabel(type);
        }

        static inline label validAddressing
        (
            addressingType type,
            label typeI
        )
        {
            return !invalidAddressing(type, typeI);
        }

        static inline label invalidMapLabel()
        {
            return -1;
        }

        static inline label invalidMap
        (
            label typeI
        )
        {
            return typeI == invalidMapLabel();
        }

        static inline label validMap
        (
            label typeI
        )
        {
            return !invalidMap(typeI);
        }

        static inline label faceAddressingToMapI
        (
            label faceAddressingI
        )
        {
            return mag(faceAddressingI) - 1;
        }

        static inline void faceAddressingToMap
        (
            const labelList& faceAddressing,
            labelList& faceMap
        )
        {
            forAll (faceAddressing, faceI)
            {
                faceMap[faceI] =
                    faceAddressingToMapI(faceAddressing[faceI]);
            }
        }

        static inline label faceMapToAddressingI
        (
            label faceMapI,
            bool faceFlippedI
        )
        {
            label faceAddressingI = faceMapI + 1;

            return (faceFlippedI) ? -faceAddressingI : faceAddressingI;
        }

        static inline void faceMapToAddressing
        (
            labelList& faceMap,
            const boolList& faceFlipped
        )
        {
            forAll (faceMap, faceI)
            {
                faceMap[faceI] =
                    faceMapToAddressingI(faceMap[faceI], faceFlipped[faceI]);
            }
        }


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        addressingTypes(const addressingTypes&);

        //- Disallow default bitwise assignment
        void operator=(const addressingTypes&);


public:

    // Constructors

        //- Null constructor
        addressingTypes() {};

    // Destructor

        ~addressingTypes() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
