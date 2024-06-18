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
    Foam::edgeBiotSavart

Description
    ...

SourceFiles
    edgeBiotSavart.C

\*---------------------------------------------------------------------------*/

#ifndef edgeBiotSavart_H
#define edgeBiotSavart_H

#include "IOdictionary.H"
#include "NamedEnum.H"
#include "fvMesh.H"
#include "volFields.H"
#include "featureEdgeMesh.H"
#include "OFstream.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "calculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class edgeBiotSavart Declaration
\*---------------------------------------------------------------------------*/

class edgeBiotSavart
:
    public IOdictionary
{
public:

    //- Inductor data
    struct InductorData
    {
        word name;
        fileName file;
        Switch reverse;
        scalar current;
        scalar phase;
    };

    //- Real/Imaginary part enums
    enum complexPart
    {
        REAL = 0,
        IMAG = 1
    };

    //- Static data

        static const NamedEnum<complexPart, 2> complexPartNames;


private:

    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- List of inductors
        List<InductorData> inductorData_;

        //- List of inductor meshes
        List<featureEdgeMesh*> inductorMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        edgeBiotSavart(const edgeBiotSavart&);

        //- Disallow default bitwise assignment
        void operator=(const edgeBiotSavart&);

        // Core Calculation

            // Magnetic vector potential

                //- Calculate A for given points
                //  Parameter part may be used to respect
                //  corresponding phase angle
                tmp<vectorField> calcA
                (
                    const pointField& points,
                    complexPart part = REAL
                ) const;

            // Magnetic field density

                //- Calculate B for given points
                tmp<vectorField> calcB
                (
                    const pointField& points,
                    complexPart part = REAL
                ) const;


public:

    // Runtime type information
    TypeName("edgeBiotSavart");

    // Constructor

        edgeBiotSavart(const fvMesh& mesh);

    // Destructor

        virtual ~edgeBiotSavart();

    // Member Functions

        // Magnetic vector potential

            //- Calculate A for internal cell centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> internalA
            (
                complexPart part = REAL
            ) const;

            //- Calculate A for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> boundaryPatchA
            (
                const label& patchI,
                complexPart part = REAL
            ) const;

            //- Calculate A for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> boundaryPatchA
            (
                const word& patchName,
                complexPart part = REAL
            ) const;

            //- Calculate A
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<volVectorField> A
            (
                complexPart part = REAL
            ) const;

            //- Calculate A based on Biot-Savart on bc and
            //  laplacian solution for internal field.
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            void A
            (
                volVectorField& vf,
                complexPart part = REAL
            ) const;

        // Magnetic field density

            //- Calculate B for internal cell centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> internalB
            (
                complexPart part = REAL
            ) const;

            //- Calculate B for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> boundaryPatchB
            (
                const label& patchI,
                complexPart part = REAL
            ) const;

            //- Calculate B for boundary patch faces centres
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<vectorField> boundaryPatchB
            (
                const word& patchName,
                complexPart part = REAL
            ) const;

            //- Calculate B
            //  Parameter part may be used to respect
            //  corresponding phase angle.
            tmp<volVectorField> B
            (
                complexPart part = REAL
            ) const;

        //- Write inductor/mesh data to vtk
        void writeVTK() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
