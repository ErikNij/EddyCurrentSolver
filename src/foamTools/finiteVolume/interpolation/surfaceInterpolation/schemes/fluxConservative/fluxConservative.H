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
    Foam::fluxConservative

Description
    Flux-conserving interpolation scheme class, which interpolates a field psi
    such that the face flux based on its gradient grad(psi), the corresponding
    diffusivity gamma and optionally an additionally given jump flux F will
    be preserved for each face according to:

        (gamma_f*grad(psi)_f) * n_f + F/|S_f| = const.

    The corresponding interpolation weights read:

        w = gamma_P / linearInterpolate(gamma) * f_x

    where f_x = fN/PN is the default (linear) mesh weight. In case of F=0,
    the face value is then simply:

        psi_f = w*psi_P + (1 - w)*psi_N

    If the jump flux F is present, a "correction" is necessary to get the
    right face value according to

        psi_f += (d_f & n_f) * (gamma_N - gamma_P)/linearInterpolate(gamma)
               * f_x*(1 - f_x) * |d_f|/|S_f| * F

    NOTE: The diffusivity gamma is regarded as scalar field which is constant
          in each cell and which may be discontinous at every face of the mesh.
          Mesh-skewness is currently not considered as this could not be
          realised correctly based on weights (skewness is not necessarily
          symmetric)!

Author
    Pascal Beckstein

SourceFiles
    fluxConservative.C

\*---------------------------------------------------------------------------*/

#ifndef fluxConservative_H
#define fluxConservative_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluxConservative Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxConservative
:
    public surfaceInterpolationScheme<Type>
{

private:

    // Private typedefs

        //- Surface fields
        typedef GeometricField<Type, fvsPatchField, surfaceMesh>
        surfaceTypeField;

        //- Volume fields
        typedef GeometricField<Type, fvPatchField, volMesh>
        volTypeField;

    // Private data

        //- Name of the (scalar) diffusivity
        const volScalarField& gamma_;

        //- Jump flux field
        const surfaceTypeField* jumpFluxPtr_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const fluxConservative&);


public:

    //- Runtime type information
    TypeName("fluxConservative");


    // Constructors

        //- Construct from mesh and Istream
        fluxConservative
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_
            (
                mesh.thisDb().lookupObject<volScalarField>
                (
                    word(is)
                )
            ),
            jumpFluxPtr_(NULL)
        {
            if (!is.eof())
            {
                const surfaceTypeField& jumpFlux =
                    mesh.thisDb().lookupObject<surfaceTypeField>
                    (
                        word(is)
                    );

                jumpFluxPtr_ = &jumpFlux;
            }

// TODO: Check for min(gamma) > 0
        }


        //- Construct from mesh, faceFlux and Istream
        fluxConservative
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_
            (
                mesh.thisDb().lookupObject<volScalarField>
                (
                    word(is)
                )
            ),
            jumpFluxPtr_(NULL)
        {
            if (!is.eof())
            {
                const surfaceTypeField& jumpFlux =
                    mesh.thisDb().lookupObject<surfaceTypeField>
                    (
                        word(is)
                    );

                jumpFluxPtr_ = &jumpFlux;
            }

// TODO: Check for min(gamma) > 0
        }


        //- Construct from mesh, gamma and jump flux pointer
        fluxConservative
        (
            const fvMesh& mesh,
            const volScalarField& gamma,
            const surfaceTypeField* jumpFluxPtr = NULL
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            gamma_(gamma),
            jumpFluxPtr_(jumpFluxPtr)
        {}


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return jumpFluxPtr_;
        }

        //- Return real deltaCoeffs (not any relaxed version!)
        tmp<surfaceScalarField> deltaCoeffs() const;

        //- Return real corrVecs (not any relaxed version!)
        tmp<surfaceVectorField> correctionVectors() const;

        //- Cosinus factors for one-sided gradients
        tmp<surfaceScalarField> cosAlpha() const;

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const volTypeField& vf
        ) const;

        //- Return the explicit correction to the face-interpolate
        virtual tmp<surfaceTypeField>
        correction
        (
            const volTypeField& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fluxConservativeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
