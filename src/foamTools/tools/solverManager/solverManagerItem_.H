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
    solverManagerItem

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerItem__H
#define solverManagerItem__H

#include "solverManagerWrapper.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_ITEM_ADD_DECLARATION_START(WRAPPER, TYPE, NAME)                   \
                                                                              \
    public: class Item_##NAME                                                 \
    :                                                                         \
        public solverManagerWrappers::WRAPPER<Manager, TYPE>                  \
    {                                                                         \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<Manager, TYPE> Wrapper;            \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<Manager, TYPE>::Mesh Mesh;


# define SM_ITEM_ADD_DECLARATION_INTERFACES(NAME)                             \
                                                                              \
        private: virtual void create() const;


# define SM_ITEM_ADD_DECLARATION_INTERFACES_REGION(NAME)                      \
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
        };                                                                    \
                                                                              \
        public: virtual ThisRegion::Storage& storage() const                  \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                           \
                                                                              \
        private: Item_##NAME(const Item_##NAME&);                             \
        private: void operator=(const Item_##NAME&);                          \
                                                                              \
        public: Item_##NAME                                                   \
        (                                                                     \
            const word& name,                                                 \
            const Manager& manager,                                           \
            Mesh& mesh                                                        \
        )                                                                     \
        :                                                                     \
            Wrapper(name, manager, mesh)                                      \
        {};                                                                   \
                                                                              \
        public: virtual ~Item_##NAME() {};


# define SM_ITEM_ADD_DECLARATION_END()                                        \
                                                                              \
    };


# define SM_ITEM_ADD_DEFINITION(NAME)                                         \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P3(                                                 \
        Item_##NAME, item_##NAME, virtual,                                    \
        word(#NAME), this->manager(), this->mesh())


# define SM_ITEM_ADD_DATAACCESS(TYPE, NAME)                                   \
                                                                              \
    public: TYPE& NAME() const                                                \
    {                                                                         \
        return item_##NAME().ref();                                           \
    };                                                                        \
                                                                              \
    public: TYPE* NAME##Ptr() const                                           \
    {                                                                         \
        return item_##NAME().ptr();                                           \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONITEM_ADD(TYPE, NAME)                                        \
                                                                              \
    SM_ITEM_ADD_DECLARATION_START(RegionItem_, TYPE, NAME)                    \
    SM_ITEM_ADD_DECLARATION_INTERFACES()                                      \
    SM_ITEM_ADD_DECLARATION_INTERFACES_REGION()                               \
    SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                                \
    SM_ITEM_ADD_DECLARATION_END()                                             \
    SM_ITEM_ADD_DEFINITION(NAME)                                              \
    SM_ITEM_ADD_DATAACCESS(TYPE, NAME)


# define SM_GLOBALITEM_ADD(TYPE, NAME)                                        \
                                                                              \
    SM_ITEM_ADD_DECLARATION_START(GlobalItem_, TYPE, NAME)                    \
    SM_ITEM_ADD_DECLARATION_INTERFACES()                                      \
    SM_ITEM_ADD_DECLARATION_CONSTRUCTION(NAME)                                \
    SM_ITEM_ADD_DECLARATION_END()                                             \
    SM_ITEM_ADD_DEFINITION(NAME)                                              \
    SM_ITEM_ADD_DATAACCESS(TYPE, NAME)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                  Class solverManagerWrappers::Item_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH, class TYPE>
class Item_
:
    public solverManagerBase::Item_<MESH, TYPE>,
    public solverManagerWrappers::Base<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private typedefs

        typedef solverManagerBase::Item_<MESH, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

        typedef typename solverManagerBase::Storage<MESH> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Item_(const Item_&);

        //- Disallow default bitwise assignment
        void operator=(const Item_&);

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

        Item_
        (
            const word& name,
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Base(name),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~Item_() {};

    // Member Functions

        //- Storage
        virtual Storage& storage() const = 0;

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
                const dictionary& parentDict = storageDict();
                word itemName = this->name();

                if (parentDict.found(itemName))
                {
                    return parentDict.subDict(itemName);
                }
                else
                {
                    return this->Wrapper::emptyDict();
                }
                return storageDict();
            };
};


/*---------------------------------------------------------------------------*\
               Class solverManagerItem::RegionItem_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class TYPE>
class RegionItem_
:
    public Item_<MANAGER, typename MANAGER::Mesh, TYPE>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;


private:

    // Private typedefs

        typedef Item_<MANAGER, Mesh, TYPE> Parent;
        typedef solverManagerBase::Item_<Mesh, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Region_
            <typename MANAGER::ManagerMesh, Mesh> Region_;
        typedef typename solverManagerBase::Storage<Mesh> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionItem_(const RegionItem_&);

        //- Disallow default bitwise assignment
        void operator=(const RegionItem_&);

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

        RegionItem_
        (
            const word& name,
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(name, manager, mesh)
        {};

    // Destructor

        virtual ~RegionItem_() {};

    // Member Functions

        //- Manager

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

        //- Storage
        virtual Storage& storage() const = 0;

};


/*---------------------------------------------------------------------------*\
              Class solverManagerItem::GlobalItem_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class TYPE>
class GlobalItem_
:
    public Item_<MANAGER, typename MANAGER::ManagerMesh, TYPE>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef Item_<MANAGER, Mesh, TYPE> Parent;
        typedef solverManagerBase::Item_<Mesh, TYPE> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef typename solverManagerBase::Storage<Mesh> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalItem_(const GlobalItem_&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalItem_&);

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

        GlobalItem_
        (
            const word& name,
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(name, manager, mesh)
        {};

    // Destructor

        virtual ~GlobalItem_() {};

    // Member Functions

        //- Control
        inline typename MANAGER::Control& control() const
        {
            return this->Wrapper::manager().control();
        };

        //- Settings
        inline typename MANAGER::Settings& settings() const
        {
            return this->Wrapper::manager().settings();
        };

        //- Storage
        virtual typename MANAGER::Storage& storage() const
        {
            return this->Wrapper::manager().storage();
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
