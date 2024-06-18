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
    solverManagerStorage

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerStorage_H
#define solverManagerStorage_H

#include "solverManagerWrapper.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_STORAGE_NEW_START(WRAPPER)                                        \
                                                                              \
    SM_WRAPPER_NEW(WRAPPER, Storage)


# define SM_STORAGE_NEW_INTERFACES()                                          \
                                                                              \
        private: virtual void create() const;


# define SM_STORAGE_NEW_INTERFACES_REGION()                                   \
                                                                              \
        public: virtual ThisRegion& region() const                            \
        {                                                                     \
            return static_cast<ThisRegion&>                                   \
            (                                                                 \
                manager().regions()[regionLabel]                              \
            );                                                                \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Control& control() const                   \
        {                                                                     \
            return region().control();                                        \
        };                                                                    \
                                                                              \
        public: inline ThisRegion::Settings& settings() const                 \
        {                                                                     \
            return region().settings();                                       \
        };


# define SM_STORAGE_NEW_PUBLIC()                                              \
                                                                              \
        public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_NEW()                                               \
                                                                              \
    SM_STORAGE_NEW_START(RegionStorage)                                       \
    SM_STORAGE_NEW_INTERFACES()                                               \
    SM_STORAGE_NEW_INTERFACES_REGION()                                        \
    SM_STORAGE_NEW_PUBLIC()


# define SM_GLOBALSTORAGE_NEW()                                               \
                                                                              \
    SM_STORAGE_NEW_START(GlobalStorage)                                       \
    SM_STORAGE_NEW_INTERFACES()                                               \
    SM_STORAGE_NEW_PUBLIC()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSTORAGE_ADD()                                               \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_REGION_REGIONSTORAGE_ADD()


# define SM_GLOBALSTORAGE_ADD()                                               \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_GLOBALSTORAGE_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
               Class solverManagerWrappers::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Storage
:
    public solverManagerBase::Storage<MESH>,
    public solverManagerWrappers::Base<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private typedefs

        typedef solverManagerBase::Storage<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Storage(const Storage&);

        //- Disallow default bitwise assignment
        void operator=(const Storage&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


public:

    // Constructor

        Storage
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~Storage() {};

    // Member Functions

        // Basic

            //- Mesh
            virtual MESH& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

        // Dictionaries

            //- Dictionary
            virtual const dictionary& dict() const
            {
                return storageDict();
            };
};


/*---------------------------------------------------------------------------*\
              Class solverManagerStorage::RegionStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class RegionStorage
:
    public Storage<MANAGER, typename MANAGER::Mesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;


private:

    // Private typedefs

        typedef Storage<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Storage<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Region_
            <typename MANAGER::ManagerMesh, Mesh> Region_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionStorage(const RegionStorage&);

        //- Disallow default bitwise assignment
        void operator=(const RegionStorage&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                const dictionary& parentDict = this->Wrapper::settingsDict();
                word regionName = this->Parent::mesh().name();

                if (parentDict.found(regionName))
                {
                    return parentDict.subDict(regionName);
                }
                else
                {
                    return this->Wrapper::emptyDict();
                }
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                const dictionary& parentDict = this->Wrapper::storageDict();
                word regionName = this->Parent::mesh().name();

                if (parentDict.found(regionName))
                {
                    return parentDict.subDict(regionName);
                }
                else
                {
                    return this->Wrapper::emptyDict();
                }
            };

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


public:

    // Constructor

        RegionStorage
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(manager, mesh)
        {};

    // Destructor

        virtual ~RegionStorage() {};

    // Member Functions

        //- Global control
        inline typename MANAGER::Control& globalControl() const
        {
            return this->Wrapper::manager().control();
        };

        // - Global settings
        inline typename MANAGER::Settings& globalSettings() const
        {
            return this->Wrapper::manager().settings();
        };

        //- Global storage
        inline typename MANAGER::Storage& globalStorage() const
        {
            return this->Wrapper::manager().storage();
        };

        //- Region
        virtual Region_& region() const = 0;
};


/*---------------------------------------------------------------------------*\
             Class solverManagerStorage::GlobalStorage Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class GlobalStorage
:
    public Storage<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef Storage<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Storage<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalStorage(const GlobalStorage&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalStorage&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const
            {
                return this->Wrapper::settingsDict();
            };

            //- Storage dictionary
            virtual const dictionary& storageDict() const
            {
                return this->Wrapper::storageDict();
            };

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


public:

    // Constructor

        GlobalStorage
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(manager, mesh)
        {};

    // Destructor

        virtual ~GlobalStorage() {};

    // Member Functions

        //- Control
        inline typename MANAGER::Control& control() const
        {
            return this->Wrapper::manager().control();
        };

        // Settings
        inline typename MANAGER::Settings& settings() const
        {
            return this->Wrapper::manager().settings();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
