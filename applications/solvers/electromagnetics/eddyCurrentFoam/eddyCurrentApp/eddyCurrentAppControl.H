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
    Foam::eddyCurrentApp::Control

Description
    ...

SourceFiles
    eddyCurrentAppControlI.H
    eddyCurrentAppControl.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef eddyCurrentAppControl_H
#define eddyCurrentAppControl_H

#include "eddyCurrentAppInclude.H"

#include "solutionControl.H"
#include "solverPerformanceData.H"

// TODO: Get rid of all preprocessor statements involving eddyCurrentAppLink_H.
//       But how? Using class hierarchies will get very complicated and needs
//       derived classes for Region (which is created by macros), Control and
//       finally the Manager itself. Maybe a switch to turn on/off the
//       link-functionality would be more appropriate...

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace eddyCurrentApp
{

/*---------------------------------------------------------------------------*\
                   Class eddyCurrentApp::Control Declaration
\*---------------------------------------------------------------------------*/

class Control
:
    public solutionControl
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        Control(const Control&);

        //- Disallow default bitwise assignment
        void operator=(const Control&);


protected:

#ifdef eddyCurrentAppLink_H

    // Protected structs

        //- Electromagnetic settings
        struct emSettings
        {
            Switch enabled;
            int outputTimeIndexCycle;
            int timeIndexCycle;
            scalar timeCycle;
            scalar relDeltaAmax;
        };

        //- Electromagnetic update data
        struct emUpdateData
        {
            Switch update;
            int counter;
            int outputTimeIndex;
            scalar lastTime;
        };

#endif

    // Protected data

        //- Mesh
        const regionFvMesh& mesh_;

        // Setup

            //- Is mesh 3D?
            const bool meshIs3D_;

            //- Solution direction
            const Vector<label> solutionDir_;

            //- Name of conductor interface
            const word interfacePatchName_;

            //- Label of conductor interface
            const label interfacePatchLabel_;

            //- Maximum number of EDDYCURRENT correctors
            label nCorrEDDYCURRENT_;

            //- Maximum number of EDDYCURRENT sub-correctors
            label nSubCorrEDDYCURRENT_;

            //- Switch
            Switch switchV_;

            //- Tolerance scale
            HashTable<scalar> tolScales_;

        // Sub-iteration

            // Evolution

                //- Sub-iteration loop index
                label subCorr_;

                //- Sub-iteration scale factor
                scalar subScale_;

#ifdef eddyCurrentAppLink_H

        // Update

            //- Electromagnetic settings dictionary
            const dictionary& emUpdateSettingsDict_;

            //- Electromagnetic settings
            emSettings emSettings_;

            //- Electromagnetic update dictionary
            IOdictionary emUpdateDataDict_;

            //- Electromagnetic update data
            emUpdateData emUpdateData_;

            //- Electromgnetic update state for relDeltaA
            bool emUpdateRelDeltaA_;

            //- Electromgnetic update init state for relDeltaA
            bool emUpdateRelDeltaAInitialized_;

#endif

    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        // Sub-iteration

            //- Return true if all sub-convergence checks are satisfied
            bool subCriteriaSatisfied();

#ifdef eddyCurrentAppLink_H

        // Update

            //- Counter is zero
            inline bool updateZeroCounter() const;

            //- Output time cycle interval is due
            inline bool updateOutputTimeIndex() const;

            //- Time index cycle interval is due
            inline bool updateTimeIndex() const;

            //- Time cycle interval is due
            inline bool updateTime() const;

            //- Relative A deviation has grown to large
            bool updateRelDeltaA(label movedRegionI = -1);

            //- Update is due
            bool update(label movedRegionI = -1);

#endif


public:

    //- Run-time type information
    TypeName("eddyCurrentApp::Control");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        Control
        (
            regionFvMesh& mesh,
            const word& dictName = "EDDYCURRENT"
        );


    //- Destructor
    virtual ~Control() {};


    // Member Functions

        // Setup

            //- Mesh geometry is 3D?
            inline bool meshIs3D() const;

            //- Mesh solution direction
            inline const Vector<label>& solutionDir() const;

            //- Reference to conductor interface name
            inline const word& interfaceName() const;

            //- Return conductor interface ID
            inline label interfaceLabel() const;

            //- Switch for rotatory solution of VRe/VIm
            inline bool switchV();

        // Evolution

            //- Main control loop
            virtual bool loop();

        // Sub-iteration

            // Access

                //- Current sub-corrector loop index
                inline label subCorr() const;

                //- Return solution dictionary of sub-iteration
                dictionary subDict
                (
                    label regionI,
                    const word& name
                ) const;

                //- Set tolerance scale
                void setTolScale
                (
                    const word& name,
                    const scalar& tolScale = 1.0
                );

            // Evolution

                //- Sub control loop
                bool subLoop();

#ifdef eddyCurrentAppLink_H

        // Update

            //- Check if update is due
            bool needsUpdate(label movedRegionI = -1);

#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eddyCurrentApp

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eddyCurrentAppControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyCurrentAppControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
