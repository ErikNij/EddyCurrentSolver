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
    Foam::fv::jumpGaussLaplacian

Description
    ...

SourceFiles
    jumpGaussLaplacian.C

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#ifndef jumpGaussLaplacian_H
#define jumpGaussLaplacian_H

#include "laplacianScheme.H"
#include "fluxConservative.H"
#include "harmonic.H"
#include "snGradScheme.H"
#include "uncorrectedSnGrad.H"

// TODO: Two-sided snGradSchemes necessary!

// TODO: Only include headers which are really needed!

// TODO: The diffusivity gamma is currently regarded as
//       discontinous at every face!

// TODO: Proper template specialisation? The diffusivity may be
//       any GType. Specialization is necessary for e.g. tensors. Currently
//       we simply use the mag-function as workaround and in order to avoid
//       having to create these.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class jumpGaussLaplacian Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class jumpGaussLaplacian
:
    public fv::laplacianScheme<Type, GType>
{

private:

    // Private typedefs

        //- Surface gamma fields
        typedef GeometricField<GType, fvsPatchField, surfaceMesh>
        surfaceGTypeField;

        //- Volume gamma fields
        typedef GeometricField<GType, fvPatchField, volMesh>
        volGTypeField;

        //- Grad type
        typedef typename outerProduct<vector, Type>::type gradType;

        //- Surface fields
        typedef GeometricField<Type, fvsPatchField, surfaceMesh>
        surfaceTypeField;

        //- Volume fields
        typedef GeometricField<Type, fvPatchField, volMesh>
        volTypeField;

        //- Volume grad fields
        typedef GeometricField<gradType, fvPatchField, volMesh>
        volGradTypeField;

    // Private data

        //- Name of the (scalar) diffusivity
        const volGTypeField& gamma_;

        //- Jump flux field
        const surfaceTypeField* jumpFluxPtr_;

        //- Flux-conservative interpolation scheme for vf
        const fluxConservative<Type> vfFluxConsScheme_;

    // Private Member Functions

        //- Check for valid schmeme settings
        void checkSchemes() const;

        //- Uncorrected Laplacian
        tmp<fvMatrix<Type> > fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const volTypeField& vf
        );

        //- Add jump flux to source of matrix
        void addJumpFlux
        (
            fvMatrix<Type>& fvm,
            const surfaceScalarField& gammaMagSf,
            const volTypeField& vf
        );

        //- Add correction for both side of face to source of matrix
        void addSnGradsCorrection
        (
            fvMatrix<Type>& fvm,
            const surfaceScalarField& gammaMagSf,
            const volTypeField& vf
        );

        //- Disallow default bitwise copy construct
        jumpGaussLaplacian(const jumpGaussLaplacian&);

        //- Disallow default bitwise assignment
        void operator=(const jumpGaussLaplacian&);


public:

    //- Runtime type information
    TypeName("jumpGauss");

    // Constructors

        //- Construct from Istream
        jumpGaussLaplacian
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            laplacianScheme<Type, GType>(mesh, is),
            gamma_
            (
                mesh.thisDb().lookupObject<volGTypeField>
                (
                    word(is)
                )
            ),
            jumpFluxPtr_(NULL),
            vfFluxConsScheme_(mesh, gamma_, jumpFluxPtr_)
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

            checkSchemes();

// TODO: Check for min(gamma) > 0
        }


        //- Construct from mesh, gamma and jump flux pointer
        jumpGaussLaplacian
        (
            const fvMesh& mesh,
            const volScalarField& gamma,
            const surfaceTypeField* jumpFluxPtr = NULL
        )
        :
            laplacianScheme<Type, GType>(mesh),
            gamma_(gamma),
            jumpFluxPtr_(jumpFluxPtr),
            vfFluxConsScheme_(mesh, gamma_, jumpFluxPtr_)
        {
            checkSchemes();

// TODO: Check for min(gamma) > 0
        }

    // Destructor

        virtual ~jumpGaussLaplacian()
        {}

    // Member Functions

        // -Laplacian method with surfaceField of gamma
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const surfaceGTypeField&,
            const volTypeField&
        );

        // -Laplacian method which interpolates gamma
        virtual tmp<fvMatrix<Type> > fvmLaplacian
        (
            const volGTypeField&,
            const volTypeField&
        );

        // -NOTE: Deactivated!
        virtual tmp<volTypeField> fvcLaplacian
        (
            const volTypeField&
        );

        // -Laplacian calculation with surfaceField of gamma
// TODO
        virtual tmp<volTypeField> fvcLaplacian
        (
            const surfaceGTypeField&,
            const volTypeField&
        );

        // -Laplacian calculation which interpolates gamma
// TODO
        virtual tmp<volTypeField> fvcLaplacian
        (
            const volGTypeField&,
            const volTypeField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "jumpGaussLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
