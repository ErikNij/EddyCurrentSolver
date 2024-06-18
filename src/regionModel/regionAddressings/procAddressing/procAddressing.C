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

#include "procAddressing.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(procAddressing, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelIOList* Foam::procAddressing::readAddressing
(
    addressingType type
) const
{
    fileName meshDir = regionName()/polyMesh::meshSubDir;

    if (regionName() == polyMesh::defaultRegion)
    {
        meshDir = polyMesh::meshSubDir;
    }

    return new labelIOList
    (
        IOobject
        (
            addressingNames[type]+word("ProcAddressing"),
            time().findInstance(meshDir, "faces"),
            meshDir,
            time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
}


Foam::labelList* Foam::procAddressing::calcFaceMap() const
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

Foam::procAddressing::procAddressing
(
    const polyMesh& mesh
)
:
    addressingTypes(),
    time_(mesh.time()),
    regionName_(mesh.name()),
    faceMapPtr_()
{}


Foam::procAddressing::procAddressing
(
    const polyMesh& mesh,
    const word& regionName
)
:
    addressingTypes(),
    time_(mesh.time()),
    regionName_(regionName),
    faceMapPtr_()
{}


Foam::procAddressing::procAddressing
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

Foam::procAddressing::~procAddressing()
{
    forAll (addressingPtrs_, type)
    {
        addressingPtrs_[type].clear();
    }

    faceMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::procAddressing::typeAddressing
(
    addressingType type
) const
{
    if (addressingPtrs_[type].empty())
    {
        if (debug)
        {
            Info<< "Foam::procAddressing::typeAddressing(...) : "
                << "Read proc "
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


const Foam::labelList& Foam::procAddressing::typeMap
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
                Info<< "Foam::procAddressing::typeMap(...) : "
                    << "Calculate proc face-map"
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
