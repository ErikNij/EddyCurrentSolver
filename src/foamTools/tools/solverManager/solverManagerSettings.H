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
    solverManagerSettings

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerSettings_H
#define solverManagerSettings_H

#include "solverManagerWrapper.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_NEW_START(WRAPPER)                                       \
                                                                              \
    SM_WRAPPER_NEW(WRAPPER, Settings)


# define SM_SETTINGS_NEW_INTERFACES()                                         \
                                                                              \
        private: virtual void read() const;


# define SM_SETTINGS_NEW_INTERFACES_REGION()                                  \
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
        public: inline ThisRegion::Storage& storage() const                   \
        {                                                                     \
            return region().storage();                                        \
        };


# define SM_SETTINGS_NEW_PUBLIC()                                             \
                                                                              \
        public:


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSETTINGS_NEW()                                              \
                                                                              \
    SM_SETTINGS_NEW_START(RegionSettings)                                     \
    SM_SETTINGS_NEW_INTERFACES()                                              \
    SM_SETTINGS_NEW_INTERFACES_REGION()                                       \
    SM_SETTINGS_NEW_PUBLIC()


# define SM_GLOBALSETTINGS_NEW()                                              \
                                                                              \
    SM_SETTINGS_NEW_START(GlobalSettings)                                     \
    SM_SETTINGS_NEW_INTERFACES()                                              \
    SM_SETTINGS_NEW_PUBLIC()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_SETTINGS_MEMBER_ADD(TYPE, NAME)                                   \
                                                                              \
        public: mutable TYPE NAME;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONSETTINGS_ADD()                                              \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_REGION_REGIONSETTINGS_ADD()


# define SM_GLOBALSETTINGS_ADD()                                              \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_GLOBALSETTINGS_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                Class solverManagerWrappers::Settings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Settings
:
    public solverManagerBase::Settings<MESH>,
    public solverManagerWrappers::Base<MANAGER, MESH>
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private typedefs

        typedef solverManagerBase::Settings<MESH> Base;
        typedef solverManagerWrappers::Base<MANAGER, MESH> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Settings(const Settings&);

        //- Disallow default bitwise assignment
        void operator=(const Settings&);

        // Basic

            //- Settings dictionary
            virtual const dictionary& settingsDict() const = 0;;

            //- Storage dictionary
            virtual const dictionary& storageDict() const = 0;

        // Evolution

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        Settings
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~Settings() {};

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
                return settingsDict();
            };
};



/*---------------------------------------------------------------------------*\
             Class solverManagerSettings::RegionSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class RegionSettings
:
    public Settings<MANAGER, typename MANAGER::Mesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;


private:

    // Private typedefs

        typedef Settings<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Settings<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Region_
            <typename MANAGER::ManagerMesh, Mesh> Region_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RegionSettings(const RegionSettings&);

        //- Disallow default bitwise assignment
        void operator=(const RegionSettings&);

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

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        RegionSettings
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(manager, mesh)
        {};

    // Destructor

        virtual ~RegionSettings() {};

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
};


/*---------------------------------------------------------------------------*\
            Class solverManagerSettings::GlobalSettings Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class GlobalSettings
:
    public Settings<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef Settings<MANAGER, Mesh> Parent;
        typedef solverManagerBase::Settings<Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GlobalSettings(const GlobalSettings&);

        //- Disallow default bitwise assignment
        void operator=(const GlobalSettings&);

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

            //- Read sequence
            virtual void read() const = 0;


public:

    // Constructor

        GlobalSettings
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Parent(manager, mesh)
        {};

    // Destructor

        virtual ~GlobalSettings() {};

    // Member Functions

        //- Control
        inline typename MANAGER::Control& control() const
        {
            return this->Wrapper::manager().control();
        };

        // Storage
        inline typename MANAGER::Storage& storage() const
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
