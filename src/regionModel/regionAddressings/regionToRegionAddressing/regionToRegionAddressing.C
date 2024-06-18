
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

#include "regionToRegionAddressing.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionToRegionAddressing, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::procAddressing& Foam::regionToRegionAddressing::getProcAddressing
(
    const word& regionName
) const
{

    if (!procAddressingPtrs_.found(regionName))
    {
        if (debug)
        {
            Info<< "Foam::regionToRegionAddressing::getProcAddressing(...) : "
                << "Create proc addressing for region " << regionName
                << endl;
        }

        procAddressingPtrs_.set
        (
            regionName,
            new procAddressing(time(), regionName)
        );
    }

    return *procAddressingPtrs_[regionName];
}

const Foam::regionAddressing& Foam::regionToRegionAddressing::getRegionAddressing
(
    const word& regionName
) const
{
    if (!regionAddressingPtrs_.found(regionName))
    {
        if (debug)
        {
            Info<< "Foam::regionToRegionAddressing::getRegionAddressing(...) : "
                << "Create region addressing for region " << regionName
                << endl;
        }

        regionAddressingPtrs_.set
        (
            regionName,
            new regionAddressing(time(), regionName)
        );
    }

    return *regionAddressingPtrs_[regionName];
}


Foam::labelList* Foam::regionToRegionAddressing::calcAddressing
(
    addressingType type,
    const word& regionName
) const
{
    if (mesh().name() == regionName)
    {
        // Empty list for same regions
        return new labelList();
    }
    else
    {
        const labelList* toMapPtr = NULL;

        if (regionName == polyMesh::defaultRegion)
        {
            if (time().processorCase())
            {
                toMapPtr = &getProcAddressing(regionName).typeMap(type);
            }
            else
            {
                // In the serial case we may just return the clean region
                // addressing from this mesh here (*)
                labelList* regionAddressingPtr = new labelList
                    (
                        getRegionAddressing(mesh().name()).typeAddressing(type)
                    );

                // Clear region addressings before return
                regionAddressingPtrs_.clear();

                return regionAddressingPtr;
            }
        }
        else
        {
            toMapPtr = &getRegionAddressing(regionName).typeMap(type);
        }

        const labelList& toMap = *toMapPtr;

        Map<label> defaultToHashMap;

        // Create hash map with all valid type labels in the default serial
        // region as key and type label of the region where we would like to
        // get our addressing to
        forAll (toMap, toTypeI)
        {
            label defaultToTypeI = toMap[toTypeI];

            if (validMap(defaultToTypeI))
            {
                defaultToHashMap.insert(defaultToTypeI, toTypeI);
            }
        }


        const labelList* fromMapPtr = NULL;
        labelList*  fromSelfMapPtr = NULL;

        if (mesh().name() == polyMesh::defaultRegion)
        {
            if (time().processorCase())
            {
                fromMapPtr = &getProcAddressing(mesh().name()).typeMap(type);
            }
            else
            {
                // In the serial case we use a simple self-map here (**)
                fromSelfMapPtr = new labelList
                    (
                        size(type, mesh()),
                        invalidMapLabel()
                    );

                labelList& fromSelfMap = *fromSelfMapPtr;

                forAll (fromSelfMap, fromTypeI)
                {
                    fromSelfMap[fromTypeI] = fromTypeI;
                }

                fromMapPtr = fromSelfMapPtr;
            }
        }
        else
        {
            fromMapPtr = &getRegionAddressing(mesh().name()).typeMap(type);
        }

        const labelList& fromMap = *fromMapPtr;


        // Create a region to region map
        labelList* mapPtr = new labelList(fromMap.size(), invalidMapLabel());

        labelList& map = *mapPtr;

        // Loop over all type labels of this region, get the the corresponding
        // label in the serial default region, and try to find the latter in
        // the hash map from the other region. If found, grab the corresponding
        // label from the other region
        forAll (map, fromTypeI)
        {
            label defaultFromTypeI = fromMap[fromTypeI];

            Map<label>::iterator iter =
                defaultToHashMap.find(defaultFromTypeI);

            if (iter != defaultToHashMap.end())
            {
                label toTypeI = iter();

                map[fromTypeI] = toTypeI;
            }
        }


        // Region to region addressing is the same as region to region mapping
        // except for faces. So lets just create a new pointer as alias to make
        // that clear. Face addressing correction will be made later
        labelList* addressingPtr = mapPtr;

        labelList& addressing = *addressingPtr;

        // Convert map into addressing for faces
        if (type == FACE)
        {
            const labelList* toAddressingPtr = NULL;

            if (regionName == polyMesh::defaultRegion)
            {
                // We do not have to check whether we have a serial case
                // here as this will never be reached then. Instead, we
                // already had returned the clean region addressing from
                // above (c.f. *)
                toAddressingPtr =
                    &getProcAddressing(regionName).typeAddressing(type);
            }
            else
            {
                toAddressingPtr =
                    &getRegionAddressing(regionName).typeAddressing(type);
            }

            const labelList& toAddressing = *toAddressingPtr;


            const labelList* fromAddressingPtr = NULL;
            labelList*  fromSelfAddressingPtr = NULL;

            if (mesh().name() == polyMesh::defaultRegion)
            {
                if (time().processorCase())
                {
                    fromAddressingPtr =
                        &getProcAddressing(mesh().name()).typeAddressing(type);
                }
                else
                {
                    // In the serial case we will now re-use the self-map from
                    // above (c.f. **) and turn it into a face self-addressing.
                    // This effectively means that we shift all labels up by 1.
                    fromSelfAddressingPtr = new labelList(*fromSelfMapPtr);

                    faceMapToAddressing
                    (
                        *fromSelfAddressingPtr,
                        boolList(fromSelfMapPtr->size(), false)
                    );

                    fromAddressingPtr = fromSelfAddressingPtr;
                }
            }
            else
            {
                fromAddressingPtr =
                    &getRegionAddressing(mesh().name()).typeAddressing(type);
            }

            const labelList& fromAddressing = *fromAddressingPtr;


            // Bare in mind here, that 'map' points to the same location
            // in memory as 'addressing' (aliasing): map == addressing
            forAll (addressing, fromFaceI)
            {
                label toFaceI = map[fromFaceI];

                // Check if face exists in other region
                if (invalidMap(toFaceI))
                {
                    // Face does not exist in the other region. Replace
                    // the invalid map label with th invalid addressing
                    // label
                    addressing[fromFaceI] = invalidAddressingLabel(type);
                }
                else
                {
                    // Face exists in the other region. Orientation is
                    // beeing checked below
                    label defaultFromFaceI = fromAddressing[fromFaceI];
                    label defaultToFaceI = toAddressing[toFaceI];

                    addressing[fromFaceI] =
                        faceMapToAddressingI
                        (
                            toFaceI,
                            defaultFromFaceI != defaultToFaceI
                        );
                }
            }

            // Delete face self-addressing pointer if it was assigned
            deleteDemandDrivenData(fromSelfAddressingPtr);
        }


        // Delete self-map pointer if it was assigned
        deleteDemandDrivenData(fromSelfMapPtr);

        // Clear all proc and region addressings
        procAddressingPtrs_.clear();
        regionAddressingPtrs_.clear();

        return addressingPtr;
    }
}


Foam::labelList* Foam::regionToRegionAddressing::calcFaceMap
(
    const word& regionName
) const
{
    const labelList& faceAddressing = typeAddressing(FACE, regionName);

    labelList* faceMapPtr = new labelList
        (
            faceAddressing.size(),
            invalidMapLabel()
        );

    faceAddressingToMap(faceAddressing, *faceMapPtr);

    return faceMapPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionToRegionAddressing::regionToRegionAddressing
(
    const polyMesh& mesh
)
:
    addressingTypes(),
    time_(mesh.time()),
    mesh_(mesh),
    procAddressingPtrs_(),
    regionAddressingPtrs_(),
    addressingPtrs_(4, HashPtrTable<labelList>()),
    addressingDictPtrs_(4, autoPtr<IOdictionary>()),
    faceMapPtr_()
{
    forAll (addressingDictPtrs_, type)
    {
        autoPtr<IOdictionary>& addressingDictPtr = addressingDictPtrs_[type];

        addressingDictPtr.set
        (
            new IOdictionary
            (
                IOobject
                (
                    addressingNames[type]+word("RegionToRegionAddressing"),
                    mesh_.facesInstance(),
                    mesh_.meshSubDir,
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                dictionary()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionToRegionAddressing::~regionToRegionAddressing()
{
    procAddressingPtrs_.clear();
    regionAddressingPtrs_.clear();

    forAll (addressingPtrs_, type)
    {
        addressingPtrs_[type].clear();
    }
    forAll (addressingDictPtrs_, type)
    {
        addressingDictPtrs_[type].clear();
    }

    faceMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::regionToRegionAddressing::typeAddressing
(
    addressingType type,
    const word& regionName
) const
{
    HashPtrTable<labelList>& ptrs = addressingPtrs_[type];
    IOdictionary& dict = addressingDictPtrs_[type]();

    if (!ptrs.found(regionName))
    {
        if (dict.found(regionName))
        {
            if (debug)
            {
                Info<< "Foam::regionToRegionAddressing::typeAddressing(...) : "
                    << "Read region-to-region "
                    << addressingNames[type] << "-addressing"
                    << " (" << mesh().name() << " -> " << regionName << ")"
                    << " from dictionary"
                    << endl;
            }

            // Read addressing from dictionary
            Istream& is = dict.lookup(regionName);

            label size = readLabel(is);

            is.readBeginList("regionToRegionAddressing");

            ptrs.set
            (
                regionName,
                new labelList(size, readLabel(is)) // index 0
            );

            token it(is); // index 1, may also be end of list

            if (it.isLabel())
            {
                labelList& addressing = *ptrs[regionName];

                addressing[1] = it.labelToken();

                for (label typeI = 2; typeI < addressing.size(); typeI++)
                {
                    addressing[typeI] = readLabel(is);
                }

                is.readEndList("regionToRegionAddressing");
            }
        }
        else
        {
            if (debug)
            {
                Info<< "Foam::regionToRegionAddressing::typeAddressing(...) : "
                    << "Calculate new region-to-region "
                    << addressingNames[type] << "-addressing"
                    << " (" << mesh().name() << " -> " << regionName << ")"
                    << endl;
            }

            // Calculate new addressing
            ptrs.set
            (
                regionName,
                calcAddressing(type, regionName)
            );

            // Write addressing to dictionary
            if (mesh().name() != regionName)
            {
                dict.set(regionName, *ptrs[regionName]);
                dict.regIOobject::write();
            }
        }
    }

    return *ptrs[regionName];
}


const Foam::labelList& Foam::regionToRegionAddressing::typeMap
(
    addressingType type,
    const word& regionName
) const
{
    if (type == FACE)
    {
        HashPtrTable<labelList>& ptrs = faceMapPtr_;

        if (!ptrs.found(regionName))
        {
            if (debug)
            {
                Info<< "Foam::regionToRegionAddressing::typeMap(...) : "
                    << "Calculate new region-to-region face-map"
                    << " (" << mesh().name() << " -> " << regionName << ")"
                    << endl;
            }

            ptrs.set
            (
                regionName,
                calcFaceMap(regionName)
            );
        }

        return *ptrs[regionName];
    }
    else
    {
        return typeAddressing(type, regionName);
    }
}


// ************************************************************************* //
