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

#include "regionAddressing.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionAddressing, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelIOList* Foam::regionAddressing::readAddressing
(
    addressingType type
) const
{
    fileName meshDir = regionName()/polyMesh::meshSubDir;

    if (regionName() == polyMesh::defaultRegion)
    {
        meshDir = polyMesh::meshSubDir;

        // Empty list for default region
        return new labelIOList
        (
            IOobject
            (
                addressingNames[type]+word("RegionAddressing"),
                time().findInstance(meshDir, "faces"),
                meshDir,
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            labelList()
        );
    }
    else
    {
        // For parallel cases, addressings need to be decomposed!

        const word facesInstance = time().findInstance(meshDir, "faces");

        IOobject typeAddressingObj
        (
            addressingNames[type]+word("RegionAddressing"),
            facesInstance,
            meshDir,
            time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        IOobject typeMapObj
        (
            addressingNames[type]+word("Map"),
            facesInstance,
            meshDir,
            time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (type == BOUNDARY) typeMapObj.rename(word("patchMap"));

        // Do we have to recover the addressing from map?
        bool fromMap = false;

        if (!typeAddressingObj.headerOk() && typeMapObj.headerOk())
        {
            typeAddressingObj = typeMapObj;

            fromMap = true;
        }

        typeAddressingObj.readOpt() = IOobject::MUST_READ;

        labelIOList* addressingPtr = new labelIOList
            (
                typeAddressingObj
            );

        // Convert facemap to face addressing if recovering from map
        if (fromMap)
        {
            typeAddressingObj.rename
            (
                addressingNames[type]+word("RegionAddressing")
            );

            if (type == FACE)
            {
                labelList& addressing = *addressingPtr;

                faceMapToAddressing
                (
                    addressing,
                    boolList(addressing.size(), false)
                );
            }
        }

        return addressingPtr;
    }
}


Foam::labelList* Foam::regionAddressing::calcFaceMap() const
{
    const labelList& faceAddressing = typeAddressing(FACE);

    labelList* faceMapPtr = new labelList
        (
            faceAddressing.size(),
            invalidMapLabel()
        );

    faceAddressingToMap(faceAddressing, *faceMapPtr);

    return faceMapPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionAddressing::regionAddressing
(
    const polyMesh& mesh
)
:
    addressingTypes(),
    time_(mesh.time()),
    regionName_(mesh.name()),
    addressingPtrs_(4, autoPtr<labelIOList>()),
    faceMapPtr_()
{}


Foam::regionAddressing::regionAddressing
(
    const polyMesh& mesh,
    const word& regionName
)
:
    addressingTypes(),
    time_(mesh.time()),
    regionName_(regionName),
    addressingPtrs_(4, autoPtr<labelIOList>()),
    faceMapPtr_()
{}


Foam::regionAddressing::regionAddressing
(
    const Time& time,
    const word& regionName
)
:
    addressingTypes(),
    time_(time),
    regionName_(regionName),
    addressingPtrs_(4, autoPtr<labelIOList>()),
    faceMapPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionAddressing::~regionAddressing()
{
    forAll (addressingPtrs_, type)
    {
        addressingPtrs_[type].clear();
    }

    faceMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::regionAddressing::typeAddressing
(
    addressingType type
) const
{
    if (addressingPtrs_[type].empty())
    {
        if (debug)
        {
            Info<< "Foam::regionAddressing::typeAddressing(...) : "
                << "Read region "
                << addressingNames[type] << "-addressing"
                << " (" << regionName() << ")"
                << endl;
        }

        addressingPtrs_[type].set
        (
            readAddressing(type)
        );
    }

    return addressingPtrs_[type]();
}


const Foam::labelList& Foam::regionAddressing::typeMap
(
    addressingType type
) const
{
    if (type == FACE)
    {
        if (faceMapPtr_.empty())
        {
            if (debug)
            {
                Info<< "Foam::regionAddressing::typeMap(...) : "
                    << "Calculate region face-map"
                    << " (" << regionName() << ")"
                    << endl;
            }

            faceMapPtr_.set
            (
                calcFaceMap()
            );
        }

        return faceMapPtr_();
    }
    else
    {
        return typeAddressing(type);
    }
}


// ************************************************************************* //
