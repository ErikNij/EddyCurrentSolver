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
    solverManager

Description
    ...

SourceFiles
    solverManager.C

\*---------------------------------------------------------------------------*/

#ifndef solverManager_H
#define solverManager_H

#include "autoPtr.H"
#include "argList.H"
#include "foamTime.H"
#include "IOdictionary.H"
#include "solutionControl.H"

#include "solverManagerBase.H"

#include "solverManagerWrapper.H"
#include "solverManagerSettings.H"
#include "solverManagerItem_.H"
#include "solverManagerStorage.H"
#include "solverManagerRegion_.H"
#include "solverManagerRegions.H"

// TODO [Low]: This is a template class. Capitalize first letter!

// TODO [HIGH]: Pre-compile specialized versions of solverManager for all
//              meshes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# define SM_MANAGERSCOPE()                                                    \
                                                                              \
    typedef Manager::ManagerMesh GlobalMesh;                                  \
    typedef Manager::Control GlobalControl;                                   \
    typedef Manager::Settings GlobalSettings;                                 \
    typedef Manager::Storage GlobalStorage;                                   \
                                                                              \
    GlobalMesh& globalMesh = manager.mesh();                                  \
    GlobalControl& globalControl = manager.control();                         \
    GlobalSettings& globalSettings = manager.settings();                      \
    GlobalStorage& globalStorage = manager.storage();                         \
                                                                              \
    if (!&globalMesh) {};                                                     \
    if (!&globalControl) {};                                                  \
    if (!&globalSettings) {};                                                 \
    if (!&globalStorage) {};


# define SM_GLOBALREGIONSCOPE(ENUM)                                           \
                                                                              \
    SM_MANAGERSCOPE();                                                        \
                                                                              \
    typedef Manager::Region_##ENUM Region_;                                   \
                                                                              \
    typedef Region_::Control Control;                                         \
    typedef Region_::Settings Settings;                                       \
    typedef Region_::Storage Storage;                                         \
                                                                              \
    Region_& region = manager.regions().region_##ENUM();                      \
    Control& control = manager.regions().region_##ENUM().control();           \
    Settings& settings = manager.regions().region_##ENUM().settings();        \
    Storage& storage = manager.regions().region_##ENUM().storage();           \
                                                                              \
    if (!&region) {};                                                         \
    if (!&control) {};                                                        \
    if (!&settings) {};                                                       \
    if (!&storage) {};


# define SM_REGIONSCOPE(ENUM)                                                 \
                                                                              \
    typedef Manager::Region_##ENUM Region_;                                   \
                                                                              \
    typedef Region_::Mesh Mesh;                                               \
    typedef Region_::Control Control;                                         \
    typedef Region_::Settings Settings;                                       \
    typedef Region_::Storage Storage;                                         \
                                                                              \
    Region_& region = manager.regions().region_##ENUM();                      \
    Mesh& mesh = region.mesh();                                               \
    Control& control = manager.regions().region_##ENUM().control();           \
    Settings& settings = manager.regions().region_##ENUM().settings();        \
    Storage& storage = manager.regions().region_##ENUM().storage();           \
                                                                              \
    if (!&region) {};                                                         \
    if (!&mesh) {};                                                           \
    if (!&control) {};                                                        \
    if (!&settings) {};                                                       \
    if (!&storage) {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solverManagerName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(solverManager);


/*---------------------------------------------------------------------------*\
                         Class solverManager Declaration
\*---------------------------------------------------------------------------*/

template <class MESH>
class solverManager
:
    public solverManagerName
{
public:

    // Public typedefs

        typedef MESH ManagerMesh;
        typedef typename solverManagerBase::Mesh<ManagerMesh>::Type Mesh;

        typedef typename solverManagerBase::Messages<ManagerMesh> Messages;
        typedef typename solverManagerBase::Control Control;
        typedef typename solverManagerBase::Settings<ManagerMesh> Settings;
        typedef typename solverManagerBase::Storage<ManagerMesh> Storage;
        typedef typename
            solverManagerBase::Regions<ManagerMesh, Mesh> Regions;

    // Public data types

        //- Parameters
        class Parameters
        {
        public:

            // Public structs

                //- Time parameters
                struct Time
                {
                    bool adjustTimeStep;
                    scalar maxCo;
                    scalar maxDeltaT;
                    scalar CoNum;
                };

            // Data

                //- Time parameters
                Time time;
        };

    // Public typedefs

        typedef typename Parameters::Time TimeParameters;


private:

    // Private Data

        //- Parameters
        mutable autoPtr<Parameters> parametersPtr_;

        //- Messages
        mutable autoPtr<Messages> messagesPtr_;

        // Basic

            //- Reference to args
            const argList& args_;

            //- Reference to time
            Time& time_;

            //- Reference to mesh
            ManagerMesh& mesh_;

            //- Master manager?
            const bool master_;

        // Dictionaries

            //- Empty dictionary
            const dictionary emptyDict_;

            //- Properties dictionary
            mutable IOdictionary propertiesDict_;

        // Evolution

            //- Pre runTime phase?
            mutable bool prePhase_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverManager(const solverManager&);

        //- Disallow default bitwise assignment
        void operator=(const solverManager&);

        // Basic

            //- Error if not master manager
            void errorIfNotMaster() const;

        // Evolution

            //- Pre runTime phase?
            inline bool& prePhase() const
            {
                return prePhase_;
            };

            //- Read parameters
            void readParameters(Parameters& par) const;

            //- Read time parameters
            void readTimeParameters(TimeParameters& timePar) const;

            //- Time step setting
            void calcDeltaT(scalar& deltaT) const;

            //- Time step application
            void applyDeltaT() const;


public:

    //- Runtime type information
    TypeName("solverManager");

    // Constructor

        solverManager
        (
            const argList& args,
            Time& time,
            ManagerMesh& mesh,
            bool master = true
        );

    // Destructor

        virtual ~solverManager() {};

    // Member Functions

        //- Parameters
        const Parameters& parameters() const;

        //- Messages
        const Messages& messages() const;

        //- Regions
        virtual Regions& regions() const = 0;

        //- Control
        virtual Control& control() const = 0;

        //- Settings
        virtual Settings& settings() const = 0;

        //- Storage
        virtual Storage& storage() const = 0;

        // Basic

            //- Main arguments
            inline const argList& args() const
            {
                return args_;
            };

            //- Time
            inline Time& time() const
            {
                return time_;
            };

            //- Mesh
            inline ManagerMesh& mesh() const
            {
                return mesh_;
            };

            //- Is this the master manager?
            inline const bool& master() const
            {
                return master_;
            };

        // Dictionaries

            //- Control dictionary
            inline const dictionary& controlDict() const
            {
                return time_.controlDict();
            }

            //- Empty dictionary
            inline const dictionary& emptyDict() const
            {
                return emptyDict_;
            };

            //- Properties dictionary
            inline const dictionary& propertiesDict() const
            {
                return propertiesDict_;
            };

            inline const dictionary& regionsDict() const
            {
                return propertiesDict_.subDict("regions");
            };

            inline const dictionary& settingsDict() const
            {
                return propertiesDict_.subDict("settings");
            };

            //- Storage dictionary
            inline const dictionary& storageDict() const
            {
                return propertiesDict_.subDict("storage");
            };


        // Evolution

            //- Set time step width
            //  A return value of true means active
            virtual bool setDeltaT(scalar& deltaT) const
            {
                return false;
            };

            //- Set time parameters
            //  A return value of true means active
            virtual bool setTimeParameters(TimeParameters& timePar) const
            {
                return false;
            };

            //- Set Courant Number
            //  A return value of true means active
            virtual bool setCoNum(scalar& CoNum) const
            {
                return false;
            };

            //- Read sequence
            virtual void read() const;

            //- Init sequence
            virtual void init() const;

            //- Next sequence
            virtual void next() const = 0;

            //- Write sequence
            virtual void write() const = 0;

            //- Finalize sequence
            virtual void finalize() const = 0;

            //- Main time once wrapper
            bool once() const;

            //- Main time selected wrapper
            bool selected(const scalarList& selectedTimes) const;

            //- Main time loop wrapper
            bool loop() const;

            //- Main time run wrapper
            bool run() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverManager.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
