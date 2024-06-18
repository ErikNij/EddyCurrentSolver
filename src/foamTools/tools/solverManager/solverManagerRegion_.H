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
    solverManagerRegion

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerRegion__H
#define solverManagerRegion__H

#include "solverManagerWrapper.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// TODO [High]: Check if region is enabled or disabled!!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_NEW_START(ENUM)                                            \
                                                                              \
    SM_WRAPPER_NEW(Region_, Region_##ENUM)                                    \
                                                                              \
        public: static const Region::Labels regionLabel =                     \
            Region::ENUM;                                                     \
                                                                              \
        public: typedef Region_##ENUM ThisRegion;                             \
                                                                              \
        class Settings;                                                       \
        class Storage;


# define SM_REGION_NEW_INTERFACES(ENUM)                                       \
                                                                              \
        public: virtual label regionIndex() const                             \
        {                                                                     \
            return regionLabel;                                               \
        };                                                                    \
                                                                              \
        public: virtual Manager::Regions& regions() const                     \
        {                                                                     \
            return manager().regions();                                       \
        };


# define SM_REGION_NEW(ENUM)                                                  \
                                                                              \
    SM_REGION_NEW_START(ENUM)                                                 \
    SM_REGION_NEW_INTERFACES(ENUM)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_REGIONCONTROL_ADD_GLOBALCONTROL()                          \
                                                                              \
        public: typedef Manager::Control Control;                             \
                                                                              \
        public: virtual Control& control() const                              \
        {                                                                     \
            return this->manager().control();                                 \
        };

# define SM_REGIONCONTROL_ADD_GLOBALCONTROL()                                 \
                                                                              \
        SM_REGION_REGIONCONTROL_ADD_GLOBALCONTROL()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


# define SM_REGION_REGIONCONTROL_ADD_P1(TYPE, P1)                             \
                                                                              \
        SM_MANAGER_GLOBALCONTROL_ADD_P1(TYPE, P1)


# define SM_REGION_REGIONCONTROL_ADD_P2(TYPE, P1,P2)                          \
                                                                              \
        SM_MANAGER_GLOBALCONTROL_ADD_P2(TYPE, P1,P2)


# define SM_REGION_REGIONCONTROL_ADD_P3(TYPE, P1,P2,P3)                       \
                                                                              \
        SM_MANAGER_GLOBALCONTROL_ADD_P3(TYPE, P1,P2,P3)


# define SM_REGION_REGIONCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)                    \
                                                                              \
        SM_MANAGER_GLOBALCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)


# define SM_REGION_REGIONCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)                 \
                                                                              \
        SM_MANAGER_GLOBALCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_REGIONSETTINGS_ADD()                                       \
                                                                              \
        SM_WRAPPER_MEMBER_ADD_P2(                                             \
            Settings, settings, virtual,                                      \
            this->manager(), this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_REGIONSTORAGE_ADD()                                        \
                                                                              \
        SM_WRAPPER_MEMBER_ADD_P2(                                             \
            Storage, storage, virtual,                                        \
            this->manager(), this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGION_END()                                                      \
                                                                              \
    SM_WRAPPER_END()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                Class solverManagerWrappers::Region_ Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class Region_
:
    public solverManagerBase::Region_
        <typename MANAGER::ManagerMesh, typename MANAGER::Mesh>,
    public solverManagerWrappers::Base<MANAGER, typename MANAGER::Mesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::Mesh Mesh;

    // Public Member Functions

        // Linking

            // Region index
            virtual label regionIndex() const = 0;


private:

    // Private typedefs

        typedef typename MANAGER::ManagerMesh SuperMesh;

        typedef solverManagerBase::Region_ <SuperMesh, Mesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

        typedef solverManagerBase::Regions <SuperMesh, Mesh> Regions;

        typedef typename solverManagerBase::Control Control;
        typedef typename solverManagerBase::Settings<Mesh> Settings;
        typedef typename solverManagerBase::Storage<Mesh> Storage;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Region_(const Region_&);

        //- Disallow default bitwise assignment
        void operator=(const Region_&);


public:

    // Constructor

        Region_
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~Region_() {};

    // Member Functions

        //- Regions
        virtual Regions& regions() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

        // Basic

            //- Mesh
            virtual Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }

            //- Name
            inline word& name() const
            {
                return mesh().name();
            };

        // Dictionaries

            //- Settings dictionary
            inline const dictionary& settingsDict() const
            {
                const dictionary& parentDict = this->Wrapper::settingsDict();
                word regionName = this->Wrapper::mesh().name();

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
            inline const dictionary& storageDict() const
            {
                const dictionary& parentDict = this->Wrapper::storageDict();
                word regionName = this->Wrapper::mesh().name();

                if (parentDict.found(regionName))
                {
                    return parentDict.subDict(regionName);
                }
                else
                {
                    return this->Wrapper::emptyDict();
                }
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
