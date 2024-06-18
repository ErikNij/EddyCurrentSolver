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
    solverManagerRegions

Description
    ...

\*---------------------------------------------------------------------------*/

#ifndef solverManagerRegions_H
#define solverManagerRegions_H

#include "solverManagerWrapper.H"

// TODO [Low]: These are all template classes. Capitalize first letters!

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_NEW_START()                                               \
                                                                              \
    SM_WRAPPER_NEW(Regions, Regions)


# define SM_REGIONS_NEW_INTERFACES()                                          \
                                                                              \
        private: virtual void create() const;


# define SM_REGIONS_NEW()                                                     \
                                                                              \
    SM_REGIONS_NEW_START()                                                    \
    SM_REGIONS_NEW_INTERFACES()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_REGION_ADD(ENUM, MESH)                                    \
                                                                              \
        SM_WRAPPER_MEMBER_ADD_P2(                                             \
            Region_##ENUM, region_##ENUM, ,                                   \
            this->manager(), MESH)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_REGIONS_ADD()                                                     \
                                                                              \
    SM_WRAPPER_END()                                                          \
    SM_MANAGER_REGIONS_ADD()


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solverManagerWrappers
{

/*---------------------------------------------------------------------------*\
                Class solverManagerWrappers::Regions Declaration
\*---------------------------------------------------------------------------*/

template <class MANAGER>
class Regions
:
    public solverManagerBase::Regions
        <typename MANAGER::ManagerMesh, typename MANAGER::Mesh>,
    public solverManagerWrappers::Base<MANAGER, typename MANAGER::ManagerMesh>
{
public:

    // Public typedefs

        typedef typename MANAGER::ManagerMesh Mesh;


private:

    // Private typedefs

        typedef typename MANAGER::Mesh SubMesh;

        typedef solverManagerBase::Regions <Mesh, SubMesh> Base;
        typedef solverManagerWrappers::Base<MANAGER, Mesh> Wrapper;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Regions(const Regions&);

        //- Disallow default bitwise assignment
        void operator=(const Regions&);

        // Evolution

            //- Create sequence
            virtual void create() const = 0;


public:

    // Constructor

        Regions
        (
            const MANAGER& manager,
            Mesh& mesh
        )
        :
            Base(),
            Wrapper(manager, mesh)
        {};

    // Destructor

        virtual ~Regions() {};

    // Member Functions

        // Basic

            //- Mesh
            virtual Mesh& mesh() const
            {
                return this->Wrapper::wrappedMesh();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solverManagerWrappers

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
