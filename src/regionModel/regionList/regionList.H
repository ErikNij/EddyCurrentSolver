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
    Foam::regionList

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef regionList_H
#define regionList_H

#include "wordList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class regionList Declaration
\*---------------------------------------------------------------------------*/

class regionList
{
private:

    // Pricate Data

        // Region names
        wordList regionNames_;

        //-Region indices
        HashTable<label> regionIndices_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionList(const regionList&);


public:

    // Constructors

        //- Null constructor
        regionList()
        :
            regionNames_(),
            regionIndices_()
        {};

        //- Construct from region names
        inline regionList
        (
            const wordList& names
        )
        :
            regionNames_(names),
            regionIndices_()
        {
            forAll (regionNames_, regionI)
            {
                regionIndices_.set
                (
                    regionNames_[regionI],
                    regionI
                );
            }
        };

    // Destructor

        ~regionList() {};

    // Member Functions

        //- Return region count
        inline int size() const
        {
            return regionNames_.size();
        };

        //- Return all region names
        inline const wordList& names() const
        {
            return regionNames_;
        };

        //- Return region name for given index
        inline const word& name(label index) const
        {
            return regionNames_[index];
        };

        //- Return hash table with indices
        inline const HashTable<label>& indices() const
        {
            return regionIndices_;
        };

        //- Return region index for given name
        inline label index(const word& name) const
        {
            return regionIndices_[name];
        };

    // Member operators

        //- Assignment
        inline void operator=(const regionList& rl)
        {
            regionNames_ = rl.names();
            regionIndices_ = rl.indices();
        };

        //- Return region name for given index
        inline const word& operator[](label index) const
        {
            return name(index);
        };

        //- Return region index for given name
        inline label operator[](const word& name) const
        {
            return index(name);
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
