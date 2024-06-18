/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "edgeBiotSavart.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Foam::edgeBiotSavart, 0);
}


template<>
const char* Foam::NamedEnum<Foam::edgeBiotSavart::complexPart, 2>::names[] =
{
    "real",
    "imaginary"
};


const Foam::NamedEnum<Foam::edgeBiotSavart::complexPart, 2>
    Foam::edgeBiotSavart::complexPartNames;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::calcA
(
    const pointField& points,
    complexPart part
) const
{
    tmp<vectorField> tA
    (
        new vectorField(points.size(), vector::zero)
    );

    vectorField& A = tA();

    // Loop over all inductors
    forAll (inductorMeshes_, inductorI)
    {
        const scalar current = inductorData_[inductorI].current;
        const bool reverse = inductorData_[inductorI].reverse;
        const scalar phase = inductorData_[inductorI].phase/180.0 * M_PI;
        scalar complex = int(1-part)*cos(phase) - int(part)*sin(phase);

        const scalar factor = 1.0E-7 * current * complex;

        const edgeMesh& eMesh = *inductorMeshes_[inductorI];

        const edgeList& edges = eMesh.edges();
        const pointField& edgePoints = eMesh.points();

        // Loop over all edges of inductor mesh
        forAll (edges, edgeI)
        {
            const edge& curEdge = edges[edgeI];

            point point1 = point::zero;
            point point2 = point::zero;

            if (reverse)
            {
                point1 = edgePoints[curEdge.end()];
                point2 = edgePoints[curEdge.start()];
            }
            else
            {
                point1 = edgePoints[curEdge.start()];
                point2 = edgePoints[curEdge.end()];
            }

            // Loop over all given points
            forAll(points, pointI)
            {
                vector rr1  = points[pointI] - point1;
                vector rr2  = points[pointI] - point2;
                vector r2r1 = point2 - point1;

                scalar nrr1  = mag(rr1);
                scalar nrr2  = mag(rr2);
                scalar nr2r1 = mag(r2r1);

                scalar b = 0.5/nr2r1
                            * (pow(nrr1,2) + pow(nr2r1,2) - pow(nrr2,2));

                scalar a = sqrt(pow(nrr1,2) - pow(b,2));

                scalar s1 = b/nrr1;
                scalar s2 = (nr2r1 - b)/nrr2;
                scalar c1 = a/nrr1;
                scalar c2 = a/nrr2;

                A[pointI] += factor * r2r1/nr2r1
                           * ( log((c2 + s2 + 1.0) / (c2 - s2 + 1.0))
                             - log((c1 - s1 + 1.0) / (c1 + s1 + 1.0)) );
            }
        }
    }

    return tA;
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::calcB
(
    const pointField& points,
    complexPart part
) const
{
    tmp<vectorField> tB
    (
        new vectorField(points.size(), vector::zero)
    );

    vectorField& B = tB();

    // Loop over all inductors
    forAll (inductorMeshes_, inductorI)
    {
        const scalar current = inductorData_[inductorI].current;
        const bool reverse = inductorData_[inductorI].reverse;
        const scalar phase = inductorData_[inductorI].phase/180.0 * M_PI;
        scalar complex = int(1-part)*cos(phase) - int(part)*sin(phase);

        const scalar factor = 1.0E-7 * current * complex;

        const edgeMesh& eMesh = *inductorMeshes_[inductorI];

        const edgeList& edges = eMesh.edges();
        const pointField& edgePoints = eMesh.points();

        // Loop over all edges of inductor mesh
        forAll (edges, edgeI)
        {
            const edge& curEdge = edges[edgeI];

            point point1 = point::zero;
            point point2 = point::zero;

            if (reverse)
            {
                point1 = edgePoints[curEdge.end()];
                point2 = edgePoints[curEdge.start()];
            }
            else
            {
                point1 = edgePoints[curEdge.start()];
                point2 = edgePoints[curEdge.end()];
            }

            // Loop over all given points
            forAll(points, pointI)
            {
                vector rr1  = points[pointI] - point1;
                vector rr2  = points[pointI] - point2;
                vector r2r1 = point2 - point1;

                scalar nrr1  = mag(rr1);
                scalar nrr2  = mag(rr2);
                scalar nr2r1 = mag(r2r1);

                scalar b = 0.5/nr2r1
                         * (pow(nrr1,2) + pow(nr2r1,2) - pow(nrr2,2));

                scalar aSquare = pow(nrr1,2) - pow(b,2);

                B[pointI] += factor / aSquare/nr2r1
                           * (b/nrr1 + (nr2r1-b)/nrr2) * (rr1 ^ r2r1);
            }
        }
    }

    return tB;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeBiotSavart::edgeBiotSavart(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "edgeBiotSavartProperties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    inductorData_(0, InductorData()),
    inductorMeshes_(0, NULL)
{
    const dictionary& inDict = this->subDict("inductors");

    forAllConstIter(dictionary, inDict, iter)
    {
        InductorData curInData = InductorData();

        curInData.name = word(iter().keyword());

        const dictionary& curInDict = inDict.subDict(iter().keyword());

        curInData.file = word(curInDict.lookup("file"));
        curInData.current = readScalar(curInDict.lookup("current"));
        curInData.reverse = curInDict.lookupOrDefault("reverse",Switch(false));
        curInData.phase = curInDict.lookupOrDefault("phase", 0.0);

        inductorData_.setSize
        (
            inductorData_.size()+1,
            curInData
        );

        inductorMeshes_.setSize
        (
            inductorMeshes_.size()+1,
            new featureEdgeMesh
            (
                IOobject
                (
                    "featureEdgeMesh"/curInData.file,
                    mesh.time().constant(),
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::edgeBiotSavart::~edgeBiotSavart()
{
    forAll (inductorMeshes_, inductorI)
    {
        deleteDemandDrivenData(inductorMeshes_[inductorI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::internalA
(
    complexPart part
) const
{
    Info << "Biot-Savart for A on internal field"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return calcA(mesh_.C().internalField(), part);
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::boundaryPatchA
(
    const label& patchI,
    complexPart part
) const
{
    Info << "Biot-Savart for A on boundary field"
        << " of patch " << "'" << mesh_.boundaryMesh()[patchI].name() << "'"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return calcA(mesh_.C().boundaryField()[patchI], part);
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::boundaryPatchA
(
    const word& patchName,
    complexPart part
) const
{
    label patchI = mesh_.boundaryMesh().findPatchID(patchName);

    Info << "Biot-Savart for A on boundary field "
        << " of patch " << "'" << patchName << "'"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return calcA(mesh_.C().boundaryField()[patchI], part);
}


Foam::tmp<Foam::volVectorField> Foam::edgeBiotSavart::A
(
    complexPart part
) const
{
    tmp<volVectorField> tA
    (
        new volVectorField
        (
            IOobject
            (
                "A",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );

    volVectorField& A = tA();

    forAll (A.boundaryField(), patchI)
    {
        const fvPatchVectorField& pvf = A.boundaryField()[patchI];

        if (isA<calculatedFvPatchVectorField>(pvf))
        {
            A.boundaryField()[patchI] = boundaryPatchA(patchI, part);
        }
    }

    A.internalField() = internalA(part);

    return tA;
}


void Foam::edgeBiotSavart::A
(
    volVectorField& vf,
    complexPart part
) const
{
    tmp<volVectorField> tA
    (
        new volVectorField
        (
            IOobject
            (
                vf.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimLength,
                vector::zero
            ),
            fixedValueFvPatchVectorField::typeName
        )
    );

    volVectorField& A = tA();

    forAll (A.boundaryField(), patchI)
    {
        const fvPatchVectorField& pvf = A.boundaryField()[patchI];

        if (isA<fixedValueFvPatchVectorField>(pvf))
        {
            A.boundaryField()[patchI] == boundaryPatchA(patchI, part);
        }
    }

    int nCorr = readInt(this->lookup("nNonOrthogonalCorrectors"));

    for (int corr=0; corr<nCorr; corr++)
    {
        fvVectorMatrix AEqn(fvm::laplacian(A));

        AEqn.solve();
    }

    vf == tA;
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::internalB
(
    complexPart part
) const
{
    Info << "Biot-Savart for B on internal field"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return calcB(mesh_.C().internalField(), part);
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::boundaryPatchB
(
    const label& patchI,
    complexPart part
) const
{
    Info << "Biot-Savart for B on boundary field"
        << " of patch " << "'" << mesh_.boundaryMesh()[patchI].name() << "'"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return calcB(mesh_.C().boundaryField()[patchI], part);
}


Foam::tmp<Foam::vectorField> Foam::edgeBiotSavart::boundaryPatchB
(
    const word& patchName,
    complexPart part
) const
{
    label patchI = mesh_.boundaryMesh().findPatchID(patchName);

    Info << "Biot-Savart for B on boundary field "
        << " of patch " << "'" << patchName << "'"
        << " (" << complexPartNames[part] << ")"
        << endl;

    return boundaryPatchB(patchI, part);
}


Foam::tmp<Foam::volVectorField> Foam::edgeBiotSavart::B
(
    complexPart part
) const
{
    tmp<volVectorField> tB
    (
        new volVectorField
        (
            IOobject
            (
                "B",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                word(),
                dimVoltage*dimTime/dimArea,
                vector::zero
            ),
            calculatedFvPatchVectorField::typeName
        )
    );

    volVectorField& B = tB();

    forAll (B.boundaryField(), patchI)
    {
        const fvPatchVectorField& pvf = B.boundaryField()[patchI];

        if (isA<calculatedFvPatchVectorField>(pvf))
        {
            B.boundaryField()[patchI] = boundaryPatchB(patchI, part);
        }
    }

    B.internalField() = internalB(part);

    return tB;
}


void Foam::edgeBiotSavart::writeVTK() const
{
    // Loop over all inductors
    forAll (inductorMeshes_, inductorI)
    {
        const edgeMesh& eMesh = *inductorMeshes_[inductorI];

        const edgeList& edges = eMesh.edges();
        const pointField& edgePoints = eMesh.points();

        pointField points(edgePoints.size() + edges.size(), vector::zero);

        forAll (edgePoints, edgePointI)
        {
            const label pointI = edgePointI;

            points[pointI] = edgePoints[pointI];
        }

        forAll (edges, edgeI)
        {
            const label pointI = edgePoints.size() + edgeI;

            const edge& curEdge = edges[edgeI];

            points[pointI] = curEdge.centre(edgePoints);
        }

        const fileName meshFileName = inductorData_[inductorI].file;
        const fileName vtkFileName = mesh_.time().constant()
                                   / "featureEdgeMesh"
                                   / meshFileName.lessExt() + ".vtk";

        OFstream vtk(vtkFileName);

        // Write header
        {
            vtk << "# vtk DataFile Version 3.1" << nl
                << vtkFileName.name() << nl
                << "ASCII" << nl
                << "DATASET UNSTRUCTURED_GRID" << nl;

            vtk << nl;
        }

        // Write points
        {
            vtk << "POINTS" << ' '
                << points.size() << ' '
                << "float" << nl;

            forAll (points, pointI)
            {
                vtk << points[pointI].x() << ' '
                    << points[pointI].y() << ' '
                    << points[pointI].z() << nl;
            }

            vtk << nl;
        }

        // Write cells and cell types
        {
            vtk << "CELLS" << ' '
                << edges.size() << ' '
                << 3*edges.size() << nl;

            forAll (edges, edgeI)
            {
                const edge& curEdge = edges[edgeI];

                vtk << 2 << ' '
                    << curEdge.start() << ' '
                    << curEdge.end() << nl;
            }

            vtk << nl;

            vtk << "CELL_TYPES" << ' '
                << edges.size() << nl;

            forAll (edges, edgeI)
            {
                vtk << 3 << nl;
            }

            vtk << nl;
        }

        // Write cell data
        {
            vtk << "CELL_DATA" << ' '
                << edges.size() << nl;

            // Phase (phi0)
            {
                vtk << "SCALARS" << ' '
                    << "phi0" << ' '
                    << "float" << nl;

                vtk << "LOOKUP_TABLE default" << nl;

                forAll (edges, edgeI)
                {
                    scalar phi0I = inductorData_[inductorI].phase;

                    vtk << phi0I << nl;
                }
            }

            // Current (I0)
            {
                vtk << "VECTORS" << ' '
                    << "I0" << ' '
                    << "float" << nl;

                forAll (edges, edgeI)
                {
                    const edge& curEdge = edges[edgeI];

                    vector I0I = curEdge.vec(edgePoints)
                               / curEdge.mag(edgePoints);

                    I0I *= inductorData_[inductorI].current / mag(I0I);

                    if (inductorData_[inductorI].reverse) I0I *= -1;

                    vtk << I0I.x() << ' '
                        << I0I.y() << ' '
                        << I0I.z() << nl;
                }
            }

            vtk << nl;
        }

        // Write point data
        {
            vtk << "POINT_DATA" << ' '
                << points.size() << nl;

            // Phase (phi0)
            {
                vtk << "SCALARS" << ' '
                    << "phi0" << ' '
                    << "float" << nl;

                vtk << "LOOKUP_TABLE default" << nl;

                forAll (points, pointsI)
                {
                    scalar phi0I = inductorData_[inductorI].phase;

                    vtk << phi0I << nl;
                }
            }

            // Current (I0)
            {
                vtk << "VECTORS" << ' '
                    << "I0" << ' '
                    << "float" << nl;

                forAll (edgePoints, edgePointI)
                {
                    vector I0I = vector::zero;

                    vtk << I0I.x() << ' '
                        << I0I.y() << ' '
                        << I0I.z() << nl;
                }

                forAll (edges, edgeI)
                {
                    const edge& curEdge = edges[edgeI];

                    vector I0I = curEdge.vec(edgePoints)
                               / curEdge.mag(edgePoints);

                    I0I *= inductorData_[inductorI].current / mag(I0I);

                    if (inductorData_[inductorI].reverse) I0I *= -1;

                    vtk << I0I.x() << ' '
                        << I0I.y() << ' '
                        << I0I.z() << nl;
                }
            }

            vtk << nl;
        }
    }
}


// ************************************************************************* //
