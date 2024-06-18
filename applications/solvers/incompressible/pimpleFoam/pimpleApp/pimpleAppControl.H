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
    Foam::pimpleApp::Control

Description
    Specialization of the pimpleControl class for pimpleFoam.

SourceFiles
    pimpleAppControlI.H
    pimpleAppControl.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef pimpleAppControl_H
#define pimpleAppControl_H

#include "pimpleAppInclude.H"

#include "pimpleControl.H"
#include "solverPerformanceData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace pimpleApp
{

/*---------------------------------------------------------------------------*\
                     Class pimpleApp::Control Declaration
\*---------------------------------------------------------------------------*/

class Control
:
    public pimpleControl
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        Control(const Control&);

        //- Disallow default bitwise assignment
        void operator=(const Control&);


protected:

    // Protected structs

        //- Residual storage
        struct residualStorage
        {
            bool stored;
            scalar residual;
        };

    // Protected data

        //- Reference cell for p
        label pRefCell_;

        //- Reference value for p
        scalar pRefValue_;

        //- List of stored residuals per field
        List<residualStorage> residualStorage_;

    // Protected Member Functions

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();


public:

    //- Run-time type information
    TypeName("pimpleApp::Control");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        Control(fvMesh& mesh, const word& dictName="PIMPLE");


    //- Destructor
    virtual ~Control();


    // Member Functions

        // Access

            //- Return reference cell for p
            inline label pRefCell() const;

            //- Return reference value for p
            inline const scalar& pRefValue() const;

            //- Set reference cell and value for p from dict
            inline void setpRefCell(const volScalarField& p);

            //- Set reference cell and value for p directly
            inline void setpRefCell(label i, const scalar& v);

        // Solution control

            //- Helper function to identify first inner iteration
            inline bool firstInnerIter() const;

        // Evolution

            //- Manually store current residuals for the named field
            void storeResiduals(const word& name);

            // Stop non-orthogonal correction after the next (final)
            // iteration, if the last number of solution iterations
            // was "0" for the named field
            void skipZeroNonOrtho(const word& name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pimpleApp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleAppControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pimpleAppControl.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
