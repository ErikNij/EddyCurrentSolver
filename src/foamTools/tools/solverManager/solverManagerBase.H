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
    solverManagerBase

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerBase_H
#define solverManagerBase_H

#include "dynamicFvMesh.H"
#include "regionDynamicFvMesh.H"

#include "solutionControl.H"
#include "autoPtr.H"
#include "HashTable.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class MESH> class solverManager;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerBase
{

template <class MESH> class Mesh;
template <class MESH> class Settings;
template <class MESH, class TYPE> class Item_;
template <class MESH> class Storage;
template <class MANAGERMESH, class MESH> class Region_;
template <class MANAGERMESH, class MESH> class Regions;

/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Mesh Declaration
\*---------------------------------------------------------------------------*/

template <typename MESH>
struct Mesh
{
    typedef typename MESH::Mesh Type;
};

template <>
struct Mesh<dynamicFvMesh>
{
    typedef dynamicFvMesh Type;
};

template <>
struct Mesh<regionDynamicFvMesh>
{
    typedef dynamicFvMesh Type;
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Control Declaration
\*---------------------------------------------------------------------------*/

typedef solutionControl Control;


/*---------------------------------------------------------------------------*\
                   Class solverManagerBase::Create Declaration
\*---------------------------------------------------------------------------*/

class Create
{
private:

    // Private Data

        // Evolution

            //- Created?
            mutable Switch created_;

    // Protected Member Functions

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


protected:

    // Protected member Functions

        // Evolution

            //- Set created
            inline void setCreated() const
            {
                created_ = true;
            };

    // Constructor

        Create() : created_(false) {};

    // Destructor

        virtual ~Create() {};


public:

    // Member Functions

        // Evolution

            //- Created?
            inline bool created() const
            {
                return created_;
            };

            //- Create sequence wrapper
            inline bool checkCreate() const
            {
                if (!created_)
                {
                    created_ = true;

                    create();

                    return true;
                }

                return false;
            };
};


/*---------------------------------------------------------------------------*\
                    Class solverManagerBase::Read Declaration
\*---------------------------------------------------------------------------*/

class Read
{
private:

    // Protected Member Functions

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


protected:

    // Constructor

        Read() {};

    // Destructor

        virtual ~Read() {};


public:

    // Member Functions

        // Evolution

            //- Read sequence wrapper
            inline bool checkRead() const
            {
                read();

                return true;
            };
};


/*---------------------------------------------------------------------------*\
                   Class solverManagerBase::Init Declaration
\*---------------------------------------------------------------------------*/

class Init
{
private:

    // Private Data

        // Evolution

            //- Initalized?
            mutable Switch initialized_;

    // Protected Member Functions

        // Evolution

            //- Init sequence
            virtual void init() const = 0;


protected:

    // Protected member Functions

        // Evolution

            //- Set initialized
            inline void setInitialized() const
            {
                initialized_ = true;
            };

    // Constructor

        Init() : initialized_(false) {};

    // Destructor

        virtual ~Init() {};


public:

    // Member Functions

        // Evolution

            //- Initialized?
            inline bool initialized() const
            {
                return initialized_;
            };

            //- Init sequence wrapper
            inline bool checkInit() const
            {
                if (!initialized_)
                {
                    initialized_ = true;

                    init();

                    return true;
                }

                return false;
            };
};


/*---------------------------------------------------------------------------*\
                Class solverManagerBase::LinkSlave Declaration
\*---------------------------------------------------------------------------*/

class LinkSlave
{
private:

    // Private Data

        // Evolution

            //- Linked?
            mutable Switch linked_;

    // Protected Member Functions

        // Linking

            //- Link slave to master
            virtual void link() const = 0;

            //- Unlink slave from master
            virtual void unlink() const = 0;


protected:

    // Constructor

        LinkSlave() : linked_(false) {};

    // Destructor

        virtual ~LinkSlave() {};

    // Member Functions

        // Linking

            //- Is slave linked?
            inline bool linked() const
            {
                return linked_;
            };

            //- Link wrapper
            inline bool checkLink() const
            {
                if (!linked_)
                {
                    linked_ = true;

                    link();
                }

                return linked_;
            };

            //- Unlink wrapper
            inline bool checkUnlink() const
            {
                if (linked_)
                {
                    linked_ = false;

                    unlink();
                }

                return linked_;
            };
};


/*---------------------------------------------------------------------------*\
                Class solverManagerBase::LinkMaster Declaration
\*---------------------------------------------------------------------------*/

template <class TYPE>
class LinkMaster
{
private:

    // Protected Member Functions

        // Linking

            //- Link slave to master
            virtual void link(TYPE* slave) const = 0;

            //- Unlink slave from master
            virtual void unlink(const TYPE* slave) const = 0;


public:

    // Constructor

        LinkMaster() {};

    // Destructor

        virtual ~LinkMaster() {};
};


/*---------------------------------------------------------------------------*\
                 Class solverManagerBase::Settings Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Settings
:
    public Read
{
private:

    // Private Member Functions

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


protected:

    // Constructor

        Settings() : Read() {};

    // Destructor

        virtual ~Settings() {};


public:

    // Member Functions

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

        // Dictionaries

            //- Settings dictionary
            virtual const dictionary& dict() const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class solverManagerBase::ItemCore Declaration
\*---------------------------------------------------------------------------*/

class ItemCore
:
    public Create,
    public Init,
    public LinkSlave
{
private:

    // Private Data

        //- Name
        word name_;

        //- Locked switch
        mutable Switch locked_;

    // Private Member Functions

        // Evolution

            //- Create sequence
            virtual void create() const = 0;

            //- Init sequence
            virtual void init() const = 0;


protected:

    // Protected Member Functions

        // Linking

            //- Link item to storage
            virtual void link() const = 0;

            //- Unlink item from storage
            virtual void unlink() const = 0;

        // Item Data

            // Locking

                //- Reference to locked switch
                Switch& lockedSwitch() const
                {
                    return locked_;
                };

                //- Is item unlocked?
                inline bool unlocked() const
                {
                    return !locked_;
                };

                //- Is item locked?
                inline bool locked() const
                {
                    return locked_;
                };

                //- Lock item
                inline void lock() const
                {
                    locked_ = true;
                };

                //- Unlock item
                inline void unlock() const
                {
                    locked_ = false;
                };

    // Constructor

        ItemCore(const word& name)
        :
            Create(),
            Init(),
            LinkSlave(),
            name_(name),
            locked_(false)
        {};

    // Destructor

        virtual ~ItemCore() {};


public:

    // Member Functions

        //- Name
        inline const word& name() const
        {
            return name_;
        };
};


/*---------------------------------------------------------------------------*\
                   Class solverManagerBase::Item_ Declaration
\*---------------------------------------------------------------------------*/

template <class MESH, class TYPE>
class Item_
:
    public ItemCore
{
private:

    // Private Data

        //- Item data
        mutable autoPtr<TYPE> dataPtr_;

    // Private Member functions

        // Linking

            //- Link item to storage
            virtual void link() const
            {
                storage().link
                (
                    const_cast<Item_<MESH, TYPE>*>(this)
                );
            };

            //- Unlink item from storage
            virtual void unlink() const
            {
                storage().unlink(this);
            };

        // Item Data

            //- Make sequence
            void make() const
            {
                if (this->unlocked())
                {
                    this->lockedSwitch() = true;

                    Info << "Create " << this->name()
                        << " for region " << mesh().name()
                        << endl;

                    this->checkCreate();

                    this->lockedSwitch() = false;
                }
                else
                {
                    FatalErrorIn
                    (
                        "solverManagerBase::Item_"
                        + this->name() + "::make()"
                    )
                        << "Pointer for item is part of a circular "
                        << "reference. Check all occurances!"
                        << abort(FatalError);
                }
            };

            //- Make sequence wrapper
            void checkMake() const
            {
                if (dataPtr_.empty())
                {
                    make();

                    if (dataPtr_.empty())
                    {
                        FatalErrorIn
                        (
                            "solverManagerBase::Item_"
                            + this->name() + "::checkMake()"
                        )
                            << "Pointer for item is still unallocated "
                            << "after creation. Make sure you have set "
                            << "the pointer via set(TYPE*) in its "
                            << "create-function!"
                            << abort(FatalError);
                    }
                }
            };

            //- Init data sequence
            void initData() const
            {
                if (enabled())
                {
                    checkMake();
                }
                else
                {
                    FatalErrorIn
                    (
                        "solverManagerBase::Item_"
                        + this->name() + "::initData()"
                    )
                        << "Item is globally disabled and cannot be used."
                            << " Check your settings!"
                            << abort(FatalError);

                }
            };

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


            //- Init sequence (empty)
            virtual void init() const
            {
                // Init had nothing to do exept for data init,
                // but data init is only done on data access
                // to deal with disabled (optional) data.
                // This function is only kept here for later use.
            };


protected:

    // Protected Member functions

        // Item Data

            //- Set auto pointer content of item
            inline void set(TYPE* pointer) const
            {
                dataPtr_.set(pointer);
            };

            //- Set auto pointer of item directly
            inline void set(autoPtr<TYPE> autoPointer) const
            {
                dataPtr_ = autoPointer;
            };

            //- Return auto pointer of item
            inline autoPtr<TYPE>& get() const
            {
                return dataPtr_;
            };

    // Constructor

        Item_(const word& name) : ItemCore(name) {};

    // Destructor

        virtual ~Item_() {};


public:

    // Member Functions

        //- Storage
        virtual Storage<MESH>& storage() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

        // Dictionaries

            //- Item dictionary
            virtual const dictionary& dict() const = 0;

        // Enabling/Disabling

            //- Item is enabled
            inline bool enabled() const
            {
                return this->linked();
            };

            //- Item is disabled
            inline bool disabled() const
            {
                return !this->linked();
            };

            //- Enable item
            inline void enable() const
            {
                checkMake();

                this->checkLink();
            };

            //- Disable item
            inline void disable() const
            {
                this->checkUnlink();

                dataPtr_.clear();
            };

            //- Set state of item directly
            inline void setState(const bool& state) const
            {
                if (state)
                {
                    enable();
                }
                else
                {
                    disable();
                }
            };

        // Item Data

            // Pointer state

                //- Item is empty
                inline bool empty() const
                {
                    return dataPtr_.empty();
                };

                //- Item is valid
                inline bool valid() const
                {
                    return dataPtr_.valid();
                };

            // Access

                //- Return reference
                inline TYPE& ref() const
                {
                    initData();

                    return dataPtr_();
                };

                //- Return pointer
                inline TYPE* ptr() const
                {
                    initData();

                    return dataPtr_.operator->();
                };

                //- Return auto pointer content of item
                inline TYPE* getPtr() const
                {
                    if (valid())
                    {
                        return dataPtr_.operator->();
                    }
                    else
                    {
                        return NULL;
                    }
                };
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Storage Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class Storage
:
    public Create,
    public Init,
    public LinkMaster<ItemCore>
{
public:

    // Friends declarations

        //- Friend with Item
        template<class M, class T> friend class Item_;


protected:

    // Protected typedefs

        typedef HashTable<ItemCore*> ItemHashTable;


private:

    // Private Data

        //- Hash table
        mutable autoPtr<ItemHashTable> itemHashTable_;

    // Private Member Functions

        //- Hash table
        inline ItemHashTable& itemHashTable() const
        {
            if (itemHashTable_.empty())
            {
                itemHashTable_.set(new ItemHashTable());
            }

            return itemHashTable_();
        };

        // Linking

            //- Link/register pointer of item
            virtual void link(ItemCore* item) const
            {
                itemHashTable().insert(item->name(), item);
            };

            //- Link/register pointer of item
            virtual void unlink(const ItemCore* item) const
            {
                itemHashTable().erase(item->name());
            };

        // Evolution

            //- Create sequence
            virtual void create() const = 0;

            //- Init all regions in hash table
            virtual void init() const
            {
                this->checkCreate();

                forAllConstIter(
                    typename ItemHashTable, itemHashTable(), iter)
                {
                    iter()->checkInit();
                }
            };


protected:

    // Constructor

        Storage() : Create(), Init(), LinkMaster<ItemCore>() {};

    // Destructor

        virtual ~Storage() {};


public:

    // Member Functions

        // Items

            //- Return current hash table
            inline const ItemHashTable& hashTable() const
            {
                return itemHashTable();
            };

            //- Find existing item by named item and return ref cast to T
            template<class T>
            inline T& item(const word& name) const
            {
                return *dynamic_cast<T*>(itemHashTable()[name]);
            };

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

        // Dictionaries

            //- Storage dictionary
            virtual const dictionary& dict() const = 0;

    // Member Operators

        // Items

            //- Find and return existing named item
            virtual ItemCore& operator[](const word& name) const
            {
                return *itemHashTable()[name];
            };
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Region Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGERMESH, class MESH>
class Region_
:
    public Read,
    public Init,
    public LinkSlave
{
public:

    // Public Member Functions

        // Linking

            // Region index
            virtual label regionIndex() const = 0;


private:

    // Private Member Functions

        // Linking

            //- Link region to regions
            virtual void link() const
            {
                regions().link
                (
                    const_cast<Region_<MANAGERMESH, MESH>*>(this)
                );
            };

            //- Unlink region from regions
            virtual void unlink() const
            {
                regions().unlink(this);
            };

        // Evolution

            //- Read sequence
            virtual void read() const
            {
                settings().checkRead();
            };

            //- Init sequence
            virtual void init() const
            {
                storage().checkInit();
            };


protected:

    // Constructor

        Region_() : Read(), Init(), LinkSlave() {};

    // Destructor

        virtual ~Region_() {};


public:

    // Member Functions

        //- Regions
        virtual Regions<MANAGERMESH, MESH>& regions() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings<MESH>& settings() const = 0;

        //- Storage
        virtual Storage<MESH>& storage() const = 0;

        // Basic

            //- Mesh
            virtual MESH& mesh() const = 0;

        // Enabling/Disabling

            //- Region is enabled
            inline bool enabled() const
            {
                return this->linked();
            };

            //- Region is disabled
            inline bool disabled() const
            {
                return !this->linked();
            };

            //- Enable region
            inline void enable() const
            {
// TODO
                this->checkLink();
            };

            //- Disable region
            inline void disable() const
            {
// TODO
                this->checkUnlink();
            };

            //- Set state of region directly
            inline void setState(const bool& state) const
            {
                if (state)
                {
                    enable();
                }
                else
                {
                    disable();
                }
            };
};


/*---------------------------------------------------------------------------*\
                  Class solverManagerBase::Regions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGERMESH, class MESH>
class Regions
:
    public Create,
    public Read,
    public Init,
    LinkMaster<Region_<MANAGERMESH, MESH> >
{
public:

    // Friends declarations

        //- Friend with Region
        template<class MM, class M> friend class Region_;


private:

    // Private typedef

        //- Region list type
        typedef HashTable<Region_<MANAGERMESH, MESH>*, label> RegionHashTable;

    // Private Data

        //- Region list
        mutable autoPtr<RegionHashTable> regionHashTablePtr_;

    // Private Member Functions

        //- Region list
        inline RegionHashTable& regionHashTable() const
        {
            if (regionHashTablePtr_.empty())
            {
                regionHashTablePtr_.set(new RegionHashTable());
            }

            return regionHashTablePtr_();
        };

        // Linking

            //- Link/register pointer of region
            virtual void link
            (
                Region_<MANAGERMESH, MESH>* region
            ) const
            {
                regionHashTable().insert(region->regionIndex(), region);
            };

            //- Unlink/deregister pointer of region
            virtual void unlink
            (
                const Region_<MANAGERMESH, MESH>* region
            ) const
            {
                regionHashTable().erase(region->regionIndex());
            };

        // Evolution

            //- Create sequence
            virtual void create() const = 0;

            //- Read all regions in hash table
            virtual void read() const
            {
                this->checkCreate();

                forAllConstIter(
                    typename RegionHashTable, regionHashTable(), iter)
                {
                    iter()->checkRead();
                }
            };

            //- Init all regions in hash table
            virtual void init() const
            {
                this->checkCreate();

                forAllConstIter(
                    typename RegionHashTable, regionHashTable(), iter)
                {
                    iter()->checkInit();
                }
            };


protected:

    // Constructor

        Regions() : Create(), Read(), Init() {};

    // Destructor

        virtual ~Regions() {};


public:

    // Member Functions

        // Regions

            //- Return current hash table
            inline const RegionHashTable& hashTable() const
            {
                return regionHashTable();
            };


            //- Find existing named region and ref cast to T
            template<class T>
            inline T& region(const label& regionI) const
            {
                return *dynamic_cast<T*>(regionHashTable()[regionI]);
            };

        // Basic

            //- Mesh
            virtual MANAGERMESH& mesh() const = 0;

    // Member Operators

        // Regions

            //- Return region by index
            virtual Region_<MANAGERMESH, MESH>& operator[]
            (
                const label& regionI
            ) const
            {
                return *regionHashTable()[regionI];
            };
};


/*---------------------------------------------------------------------------*\
                 Class solverManagerBase::Massages Declaration
\*---------------------------------------------------------------------------*/

//- Messages
template <class MESH>
class Messages
{
public:

    // Friends declarations

        //- Friend with solverManager
        template <class M> friend class solverManager;


private:

    // Private Data

        //- Reference to args
        const argList& args_;

        //- Reference to time
        const Time& time_;

        //- Reference to mesh
        const MESH& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Messages(const Messages&);

        //- Disallow default bitwise assignment
        void operator=(const Messages&);


public:

    // Constructor

        Messages
        (
            const argList& args,
            const Time& time,
            const MESH& mesh
        )
        :
            args_(args), time_(time), mesh_(mesh)
        {};

    // Destructor

        ~Messages() {};

    // Member Functions

        //- Info about starting time loop
        inline void newLine() const
        {
            Info << nl;
        }

        //- Info about starting time loop
        inline void startTimeLoop() const
        {
            Info << "Starting time loop"
                << endl;
        }

        //- Info about current time step
        inline void timeStep() const
        {
            Info << "Time step = "
                << time_.timeIndex()
                << endl;
        };

        //- Print current time step width
        inline void timeDeltaT() const
        {
            Info << "deltaT = "
                <<  time_.deltaT().value()
                << endl;
        };

        //- Print current time
        inline void timeIs() const
        {
            Info << "Time = "
                << time_.value()
                << endl;
        };

        //- Print execution time
        inline void executionTime() const
        {
            Info << "ExecutionTime = "
                << scalar(time_.elapsedCpuTime()) << " s"
                << endl;
        };

        //- Print end statement
        inline void end() const
        {
            Info << "End" << nl << endl;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerBase

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
