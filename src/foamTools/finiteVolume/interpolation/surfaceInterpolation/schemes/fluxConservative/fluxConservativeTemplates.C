
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

\*---------------------------------------------------------------------------*/

#include "fluxConservative.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::fluxConservative<Type>::deltaCoeffs() const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tdeltaCoeffs
    (
        new surfaceScalarField
        (
            IOobject
            (
                "deltaCoeffs",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(word(), dimless/dimLength, 0.0)
        )
    );
    surfaceScalarField& deltaCoeffs = tdeltaCoeffs();
    scalarField& deltaCoeffsIn = deltaCoeffs.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const volVectorField& C = mesh.C();
    const vectorField& CIn = C.internalField();

    // Calculate internal deltaCoeffs
    forAll (owner, faceI)
    {
        // Cell labels
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector delta = CIn[nei] - CIn[own];

        deltaCoeffsIn[faceI] = 1.0/mag(delta);
    }

    // Calculate boundary deltaCoeffs
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        const unallocLabelList& faceCells = patch.patch().faceCells();

        const vectorField& CPatch = C.boundaryField()[patchI];

        scalarField& deltaCoeffsPatch = deltaCoeffs.boundaryField()[patchI];

        forAll (patch, faceI)
        {
            const label own = faceCells[faceI];

            vector deltaPatch = CPatch[faceI] - CIn[own];

            deltaCoeffsPatch[faceI] = 1.0/mag(deltaPatch);
        }
    }

    return tdeltaCoeffs;
}

template<class Type>
Foam::tmp<Foam::surfaceVectorField>
Foam::fluxConservative<Type>::correctionVectors() const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceVectorField> tcorrVecs
    (
        new surfaceVectorField
        (
            IOobject
            (
                "corrVecsOneSided",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(word(), dimless, vector::zero)
        )
    );
    surfaceVectorField& corrVecs = tcorrVecs();
    vectorField& corrVecsIn = corrVecs.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const volVectorField& C = mesh.C();
    const vectorField& CIn = C.internalField();
    const surfaceVectorField& Sf = mesh.Sf();
    const vectorField& SfIn = Sf.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Real delta coefficients
    tmp<surfaceScalarField> tdeltaCoeffs = this->deltaCoeffs();
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();
    const scalarField& deltaCoeffsIn = deltaCoeffs.internalField();

    // Calculate internal corrVecs
    forAll (owner, faceI)
    {
        // Cell labels
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector NfInI = SfIn[faceI]/magSfIn[faceI];

        vector DfInI = (CIn[nei] - CIn[own])
                     * deltaCoeffsIn[faceI];

        corrVecsIn[faceI] = pos(NfInI & DfInI) * (NfInI - DfInI);
    }

    // Calculate boundary corrVecs
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (patch.coupled())
        {
            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& deltaCoeffsPatch =
                deltaCoeffs.boundaryField()[patchI];
            const vectorField& CPatch = C.boundaryField()[patchI];
            const scalarField& magSfPatch = magSf.boundaryField()[patchI];
            const vectorField& SfPatch = Sf.boundaryField()[patchI];

            vectorField& corrVecsPatch = corrVecs.boundaryField()[patchI];

            forAll (patch, faceI)
            {
                const label own = faceCells[faceI];

                vector NfPatchI = SfPatch[faceI]/magSfPatch[faceI];

                vector DfPatchI = (CPatch[faceI] - CIn[own])
                                * deltaCoeffsPatch[faceI];

                corrVecsPatch[faceI] = NfPatchI - DfPatchI;
            }
        }
    }

    return tcorrVecs;
}

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::fluxConservative<Type>::cosAlpha() const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tcosAlpha
    (
        new surfaceScalarField
        (
            IOobject
            (
                "cosAlpha",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(word(), dimless, 1.0)
        )
    );
    surfaceScalarField& cosAlpha = tcosAlpha();
    scalarField& cosAlphaIn = cosAlpha.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const volVectorField& C = mesh.C();
    const vectorField& CIn = C.internalField();
    const surfaceVectorField& Sf = mesh.Sf();
    const vectorField& SfIn = Sf.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Real delta coefficients
    tmp<surfaceScalarField> tdeltaCoeffs = this->deltaCoeffs();
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();
    const scalarField& deltaCoeffsIn = deltaCoeffs.internalField();

    // Calculate internal cosAlpha
    forAll (owner, faceI)
    {
        // Cell labels
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector NfInI = SfIn[faceI]/magSfIn[faceI];

        vector DfInI = (CIn[nei] - CIn[own])
                     * deltaCoeffsIn[faceI];

        cosAlphaIn[faceI] = NfInI & DfInI;
    }

    // Calculate boundary cosAlpha
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        const unallocLabelList& faceCells = patch.patch().faceCells();

        const scalarField& deltaCoeffsPatch =
            deltaCoeffs.boundaryField()[patchI];
        const vectorField& CPatch = C.boundaryField()[patchI];
        const scalarField& magSfPatch = magSf.boundaryField()[patchI];
        const vectorField& SfPatch = Sf.boundaryField()[patchI];

        scalarField& cosAlphaPatch = cosAlpha.boundaryField()[patchI];

        forAll (patch, faceI)
        {
            const label own = faceCells[faceI];

            vector NfPatchI = SfPatch[faceI]/magSfPatch[faceI];

            vector DfPatchI = (CPatch[faceI] - CIn[own])
                            * deltaCoeffsPatch[faceI];

            cosAlphaPatch[faceI] = NfPatchI & DfPatchI;
        }
    }

    return tcosAlpha;
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::fluxConservative<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                "fluxConservativeWeightingFactors" + vf.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh ,
            dimless
        )
    );

    surfaceScalarField& w = tw();
    scalarField& wIn = w.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();

    // Gamma
    const volScalarField& gamma = gamma_;
    const scalarField& gammaIn = gamma.internalField();

    // Linear (!!!) interpolated gamma
    tmp<surfaceScalarField> tgammaf(linearInterpolate(gamma));
    const surfaceScalarField& gammaf = tgammaf();
    const scalarField& gammafIn = gammaf.internalField();

    // Calculate internal w
    forAll (owner, faceI)
    {
        label own = owner[faceI];

        wIn[faceI] = gammaIn[own]/gammafIn[faceI] * weightsIn[faceI];
    }

    // Calculate boundary w
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        scalarField& wPatch = w.boundaryField()[patchI];

        // Coupled patches
        if (patch.coupled())
        {
            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& weightsPatch = weights.boundaryField()[patchI];

            const scalarField& gammafPatch = gammaf.boundaryField()[patchI];

            forAll (patch, faceI)
            {
                const label own = faceCells[faceI];

                wPatch[faceI] =
                    gammaIn[own]/gammafPatch[faceI] * weightsPatch[faceI];
            }
        }
        else
        {
            // Boundary w for uncoupled patches are 1
            wPatch = 1;
        }
    }

    return tw;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fluxConservative<Type>::correction
(
    const volTypeField& vf
) const
{
    if (!jumpFluxPtr_)
    {
        FatalErrorIn
        (
            "tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > fvmGrad\n"
            "Foam::fluxConservative<Type>::correction\n"
            "(\n"
            "    GeometricField<Type, fvPatchField, volMesh>&"
            ")\n"
        )   << "Jump flux pointer not assigned.\n"
            << "Maybe this scheme was constructed from IStream without\n"
            << " getting a name for the jump flux field?"
            << " If this is true, ::correction(vf) must not be used"
            << " and ::corrected() should actually return 'false'!"
            << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    tmp<surfaceTypeField> tcorr
    (
        new surfaceTypeField
        (
            IOobject
            (
                "fluxConservativeCorrection" + vf.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<Type>(word(), vf.dimensions(), pTraits<Type>::zero)
        )
    );

    surfaceTypeField& corr = tcorr();
    Field<Type>& corrIn = corr.internalField();

    // Mesh addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Mesh and basic surface interpolation data
    const surfaceScalarField& weights = mesh.weights();
    const scalarField& weightsIn = weights.internalField();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& magSfIn = magSf.internalField();

    // Real delta coefficients
    tmp<surfaceScalarField> tdeltaCoeffs = this->deltaCoeffs();
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();
    const scalarField& deltaCoeffsIn = deltaCoeffs.internalField();

    // Cosinus factors for one-sided gradient
    tmp<surfaceScalarField> tcosAlpha = this->cosAlpha();
    const surfaceScalarField& cosAlpha = tcosAlpha();
    const scalarField& cosAlphaIn = cosAlpha.internalField();

    // Gamma
    const volScalarField& gamma = gamma_;
    const scalarField& gammaIn = gamma.internalField();

    // Linear (!!!) interpolated gamma
    tmp<surfaceScalarField> tgammaf(linearInterpolate(gamma));
    const surfaceScalarField& gammaf = tgammaf();
    const scalarField& gammafIn = gammaf.internalField();

    // Jump flux
    const surfaceTypeField& jumpFlux = *jumpFluxPtr_;
    const Field<Type>& jumpFluxIn = jumpFlux.internalField();

    // Calculate internal correction
    forAll (owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        scalar gammaRelDiff = (gammaIn[nei] - gammaIn[own])/gammafIn[faceI];

        scalar wPwN = weightsIn[faceI] * (1 - weightsIn[faceI]);

        scalar dByMagSf = 1.0/deltaCoeffsIn[faceI]/magSfIn[faceI];

        corrIn[faceI] = cosAlphaIn[faceI]
                      * gammaRelDiff * wPwN * dByMagSf
                      * jumpFluxIn[faceI];
    }

    // Calculate boundary correction
    // (Correction will be zero if grad(gamma)*n == 0)
    forAll (mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        // Coupled patches
        if (patch.coupled())
        {
            Field<Type>& corrPatch = corr.boundaryField()[patchI];

            const unallocLabelList& faceCells = patch.patch().faceCells();

            const scalarField& weightsPatch = weights.boundaryField()[patchI];
            const scalarField& deltaCoeffsPatch = deltaCoeffs.boundaryField()[patchI];
            const scalarField& magSfPatch = magSf.boundaryField()[patchI];

            const scalarField& cosAlphaPatch = cosAlpha.boundaryField()[patchI];

            const scalarField& gammaPatch = gamma.boundaryField()[patchI];
            const scalarField& gammafPatch = gammaf.boundaryField()[patchI];
            const Field<Type>& jumpFluxPatch = jumpFlux.boundaryField()[patchI];

            forAll (patch, faceI)
            {
                const label own = faceCells[faceI];

                scalar gammaRelDiff =
                    (gammaPatch[faceI] - gammaIn[own])/gammafPatch[faceI];

                scalar wPwN = weightsPatch[faceI] * (1 - weightsPatch[faceI]);

                scalar dByMagSf = 1.0/deltaCoeffsPatch[faceI]/magSfPatch[faceI];

                corrPatch[faceI] = cosAlphaPatch[faceI]
                                 * gammaRelDiff * wPwN * dByMagSf
                                 * jumpFluxPatch[faceI];
            }
        }
    }

    return tcorr;
}


// ************************************************************************* //
