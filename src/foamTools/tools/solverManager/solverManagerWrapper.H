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
    solverManagerWrapper

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerWrapper_H
#define solverManagerWrapper_H

#include "argList.H"
#include "foamTime.H"
#include "dictionary.H"

#include "solverManagerBase.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_NEW_START(WRAPPER, NAME)                                  \
                                                                              \
    public: class NAME                                                        \
    :                                                                         \
        public solverManagerWrappers::WRAPPER<Manager>                        \
    {                                                                         \
        public: typedef typename                                              \
            solverManagerWrappers::WRAPPER<Manager> Wrapper;                  \
        public: typedef typename Wrapper::Mesh Mesh;


# define SM_WRAPPER_NEW_CONSTRUCTOR(NAME)                                     \
                                                                              \
        private: NAME(const NAME&);                                           \
        private: void operator=(const NAME&);                                 \
                                                                              \
        public: NAME                                                          \
        (                                                                     \
            const Manager& manager,                                           \
            Mesh& mesh                                                        \
        )                                                                     \
        :                                                                     \
            Wrapper(manager, mesh)                                            \
        {};                                                                   \
                                                                              \
        public: virtual ~NAME() {};


# define SM_WRAPPER_NEW(WRAPPER, NAME)                                        \
                                                                              \
    SM_WRAPPER_NEW_START(WRAPPER, NAME)                                       \
    SM_WRAPPER_NEW_CONSTRUCTOR(NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_END()                                                     \
                                                                              \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_SET_P1(TYPE, NAME, P1)                             \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1));

# define SM_WRAPPER_MEMBER_SET_P2(TYPE, NAME, P1,P2)                          \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2));

# define SM_WRAPPER_MEMBER_SET_P3(TYPE, NAME, P1,P2,P3)                       \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3));

# define SM_WRAPPER_MEMBER_SET_P4(TYPE, NAME, P1,P2,P3,P4)                    \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4));

# define SM_WRAPPER_MEMBER_SET_P5(TYPE, NAME, P1,P2,P3,P4,P5)                 \
                                                                              \
                NAME##Ptr_.set(new TYPE(P1,P2,P3,P4,P5));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_WRAPPER_MEMBER_ADD_P1(TYPE, NAME, MOD, P1)                        \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P1(TYPE, NAME, P1)                      \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_ADD_P2(TYPE, NAME, MOD, P1,P2)                     \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P2(TYPE, NAME, P1,P2)                   \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_ADD_P3(TYPE, NAME, MOD, P1,P2,P3)                  \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P3(TYPE, NAME, P1,P2,P3)                \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_ADD_P4(TYPE, NAME, MOD, P1,P2,P3,P4)               \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P4(TYPE, NAME, P1,P2,P3,P4)             \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };

# define SM_WRAPPER_MEMBER_ADD_P5(TYPE, NAME, MOD, P1,P2,P3,P4,P5) \
                                                                              \
        private: mutable autoPtr<TYPE> NAME##Ptr_;                            \
                                                                              \
        public: MOD TYPE& NAME() const                                        \
        {                                                                     \
            if (NAME##Ptr_.empty())                                           \
            {                                                                 \
                SM_WRAPPER_MEMBER_SET_P5(TYPE, NAME, P1,P2,P3,P4,P5)          \
            }                                                                 \
                                                                              \
            return NAME##Ptr_();                                              \
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_MANAGER_GLOBALCONTROL_ADD_P1(TYPE, P1)                            \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P1(                                                 \
        Control, control, virtual, P1)

# define SM_MANAGER_GLOBALCONTROL_ADD_P2(TYPE, P1,P2)                         \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P2(                                                 \
        Control, control, virtual, P1,P2)

# define SM_MANAGER_GLOBALCONTROL_ADD_P3(TYPE, P1,P2,P3)                      \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P3(                                                 \
        Control, control, virtual, P1,P2,P3)

# define SM_MANAGER_GLOBALCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)                   \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P4(                                                 \
        Control, control, virtual, P1,P2,P3,P4)

# define SM_MANAGER_GLOBALCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)                \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P5(                                                 \
        Control, control, virtual, P1,P2,P3,P4,P5)


# define SM_GLOBALCONTROL_ADD_P1(TYPE, P1)                                    \
                                                                              \
    SM_MANAGER_GLOBALCONTROL_ADD_P1(TYPE, P1)


# define SM_GLOBALCONTROL_ADD_P2(TYPE, P1,P2)                                 \
                                                                              \
    SM_MANAGER_GLOBALCONTROL_ADD_P2(TYPE, P1,P2)


# define SM_GLOBALCONTROL_ADD_P3(TYPE, P1,P2,P3)                              \
                                                                              \
    SM_MANAGER_GLOBALCONTROL_ADD_P3(TYPE, P1,P2,P3)


# define SM_GLOBALCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)                           \
                                                                              \
    SM_MANAGER_GLOBALCONTROL_ADD_P4(TYPE, P1,P2,P3,P4)


# define SM_GLOBALCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)                        \
                                                                              \
    SM_MANAGER_GLOBALCONTROL_ADD_P5(TYPE, P1,P2,P3,P4,P5)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_MANAGER_GLOBALSETTINGS_ADD()                                      \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P2(                                                 \
        Settings, settings, virtual,                                          \
        *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_MANAGER_GLOBALSTORAGE_ADD()                                       \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P2(                                                 \
        Storage, storage, virtual,                                            \
        *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_MANAGER_REGIONS_ADD()                                             \
                                                                              \
    SM_WRAPPER_MEMBER_ADD_P2(                                                 \
        Regions, regions, virtual,                                            \
        *this, this->mesh())


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                   Class solverManagerWrappers::Base Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER, class MESH>
class Base
{
public:

    // Public typedefs

        typedef MESH Mesh;


private:

    // Private Data

        //- Reference to manager
        const MANAGER& manager_;

        // Basic

            //- Reference to mesh
            Mesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Base(const Base&);

        //- Disallow default bitwise assignment
        void operator=(const Base&);


protected:

    // Protected Member Functions

        // Basic

            //- Mesh
            inline Mesh& wrappedMesh() const
            {
                return mesh_;
            }


public:

    // Constructor

        Base
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            manager_(manager),
            mesh_(mesh)
        {};

    // Destructor

        virtual ~Base() {};

    // Member Functions

        //- Manager
        inline const MANAGER& manager() const
        {
            return manager_;
        };

        // Basic

            //- Main arguments
            inline const argList& args() const
            {
                return manager().args();
            };

            //- Time
            inline Time& time() const
            {
                return manager().time();
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return manager().master();
            };

        // Dictionaries

            //- Control dictionary
            inline const dictionary& controlDict() const
            {
                return time().controlDict();
            }

            //- Empty dictionary
            inline const dictionary& emptyDict() const
            {
                return manager().emptyDict();
            };

            //- Properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return manager().propertiesDict();
            };

            inline const dictionary& regionsDict() const
            {
                return manager().regionsDict();
            };

            //- Settings dictionary
            inline const dictionary& settingsDict() const
            {
                return manager().settingsDict();
            };

            //- Storage dictionary
            inline const dictionary& storageDict() const
            {
                return manager().storageDict();
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
