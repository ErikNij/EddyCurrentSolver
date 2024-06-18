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

#include "trackedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void trackedSurface::writeVTK() const
{
// TODO: Region name! Use DB().timePath()/region_name if
//       a region other than region0 is used!

    aMesh().patch().writeVTK
    (
        DB().timePath()/prefix_,
        aMesh().patch(),
        aMesh().patch().points()
    );
}


void trackedSurface::writeVTKpoints
(
    const word fieldName,
    const vectorField& pf
) const
{
    word FieldName =
        word(toupper(fieldName[0]))
      + word(fieldName.substr(1));

    // Write patch and points into VTK
    fileName name(DB().timePath()/prefix_+FieldName);
    OFstream mps(name + ".vtk");

    mps << "# vtk DataFile Version 2.0" << nl
        << name << ".vtk" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << pf.size() << " float" << nl;

    forAll (pf, pointI)
    {
        mps << pf[pointI].x() << ' '
            << pf[pointI].y() << ' '
            << pf[pointI].z() << nl;
    }

    // Write vertices
    mps << "VERTICES " << pf.size() << ' '
        << pf.size()*2 << nl;

    forAll (pf, pointI)
    {
        mps << 1 << ' ' << pointI << nl;
    }
}


void trackedSurface::writeVTKControlPoints()
{
    writeVTKpoints
    (
        word("ControlPoints"),
        controlPoints()
    );
}


template<class Type>
void
trackedSurface::writeVol
(
    const GeometricField<Type, faPatchField, areaMesh>& af
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "vol" + word(toupper(af.name()[0]))
                    + word(af.name().substr(1)),
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensioned<Type>
            (
                "0",
                af.dimensions(),
                pTraits<Type>::zero
            ),
            fixedValueFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    forAll (mesh().boundaryMesh(), patchI)
    {
        vf.boundaryField()[patchI] == pTraits<Type>::zero;
    }

    vf.boundaryField()[aPatchID()] == af;

    vf.write();

    tvf.clear();
}


void trackedSurface::writeA()
{
    // faceAreaNormals
    aMesh().faceAreaNormals().write();

    // faceCurvatures
    curvature().write();

    // Us
    Us().write();

    // fac::div(Us())
    fac::div(Us())().write();

    // fac::grad(Us())
    fac::grad(Us())().write();

    // faceCurvaturesDivNormals
    (fac::div(aMesh().faceAreaNormals()))().write();
}


void trackedSurface::writeVolA()
{
    // Us
    writeVol(Us());

    // fac::div(Us())
    writeVol((fac::div(Us()))());

    // fac::grad(Us())
    writeVol((fac::grad(Us()))());

    // faceAreaNormals
    writeVol(aMesh().faceAreaNormals());

    // faceCurvatures
    writeVol(curvature());

    // faceCurvaturesSmoothed
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed",
            fac::average(curvature())
        )
    );

    // faceCurvaturesSmoothed2
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed2",
            fac::average(fac::average(curvature()))
        )
    );

    // faceCurvaturesSmoothed3
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed3",
            fac::average(fac::average(fac::average(curvature())))
        )
    );

    // faceCurvaturesSmoothed4
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed5",
            fac::average(fac::average(fac::average(fac::average(curvature()))))
        )
    );

    // faceCurvaturesSmoothed5
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed5",
            fac::average(fac::average(fac::average(fac::average(fac::average(curvature())))))
        )
    );

    // faceCurvaturesSmoothed10
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesSmoothed10",
            fac::average(fac::average(fac::average(fac::average(fac::average(fac::average(fac::average(fac::average(fac::average(fac::average(curvature()))))))))))
        )
    );

// TODO: WTF? Why is this not working as expected?
    // faceCurvaturesFromDivOfNormals
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesFromDivOfNormals",
            -fac::div(aMesh().faceAreaNormals())
        )
    );

    // faceCurvaturesFromSPGradOfNormals
    writeVol
    (
        areaScalarField
        (
            "faceCurvaturesFromSPGradOfNormals",
            -tr(fac::grad(aMesh().faceAreaNormals()))
        )
    );

    // surfaceTension
    writeVol(surfaceTension());

    // surfaceTensionCurvatureN
    writeVol(
        areaVectorField(
            "surfaceTensionCurvatureN",
            surfaceTension()*curvature()*aMesh().faceAreaNormals()
        )
    );

    // surfaceTension
    writeVol(
        areaVectorField(
            "surfaceTensionGrad",
            fac::grad(surfaceTension())
        )
    );

    // surfaceTensionForce
    writeVol(surfaceTensionForce());

    const vectorField& nA = aMesh().faceAreaNormals().internalField();

    // normal surfaceTensionForce
    areaVectorField surfaceTensionForceNorm(
        surfaceTensionForce().name() + "Norm",
        surfaceTensionForce()
    );
    surfaceTensionForceNorm.internalField() = ((nA*nA)&surfaceTensionForce());
    writeVol(surfaceTensionForceNorm);

    // tangential surfaceTensionForce
    areaVectorField surfaceTensionForceTang(
        surfaceTensionForce().name() + "Tang",
        surfaceTensionForce()
    );
    surfaceTensionForceTang.internalField() = ((I-nA*nA)&surfaceTensionForce());
    writeVol(surfaceTensionForceTang);

    // H
    tmp<areaVectorField> tH
    (
        new areaVectorField
        (
            IOobject
            (
                "H",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedVector
            (
                word(),
                dimless,
                vector::zero
            ),
            zeroGradientFaPatchVectorField::typeName
        )
    );
    areaVectorField& H = tH();

    H.internalField() =
        controlPoints() - aMesh().areaCentres().internalField();
    H.correctBoundaryConditions();

    writeVol(H);

    // H
    tmp<areaScalarField> tmagH
    (
        new areaScalarField
        (
            IOobject
            (
                "magH",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& magH = tmagH();

    magH.internalField() =
        mag(H.internalField())
      * sign(H.internalField()&facesDisplacementDir());
    magH.correctBoundaryConditions();

    writeVol(magH);
    tmagH.clear();
    tH.clear();

    // nGradUn
    tmp<areaScalarField> tnGradUn
    (
        new areaScalarField
        (
            IOobject
            (
                "nGradUn",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& anGradUn = tnGradUn();

    anGradUn.internalField() = nGradUn();
    anGradUn.correctBoundaryConditions();

    writeVol(anGradUn);
    tnGradUn.clear();

    // phi
    tmp<areaScalarField> taPhi
    (
        new areaScalarField
        (
            IOobject
            (
                "aPhi",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& aPhi = taPhi();

    aPhi.internalField() =
        phi().boundaryField()[aPatchID()];
    aPhi.correctBoundaryConditions();

    writeVol(aPhi);

    // meshPhi
    tmp<areaScalarField> taPhiMesh
    (
        new areaScalarField
        (
            IOobject
            (
                "aPhiMesh",
                DB().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            aMesh(),
            dimensionedScalar
            (
                word(),
                dimless,
                0
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& aPhiMesh = taPhiMesh();

    aPhiMesh.internalField() =
        fvc::meshPhi(rho(),U())().boundaryField()[aPatchID()];
    aPhiMesh.correctBoundaryConditions();

    writeVol(aPhiMesh);

    // diffPhi
    writeVol
    (
        areaScalarField
        (
            "aPhiDiff",
            aPhi-aPhiMesh
        )
    );
    taPhi.clear();
    taPhiMesh.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
