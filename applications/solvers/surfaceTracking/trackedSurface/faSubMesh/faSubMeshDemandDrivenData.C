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

#include "faSubMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faSubMesh::makeFaceSubToBaseAreaMap() const
{
    if (debug)
    {
        Info<< "faSubMesh::makeFaceCurvatures() : "
            << "making sub to base face map for finite area meshes"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (faceSubToBaseAreaMapPtr_)
    {
        FatalErrorIn("faSubMesh::makeFaceCurvatures()")
            << "sub to base face map field already exists"
                << abort(FatalError);
    }

    const labelList& baseAreaFaceLabels = baseAreaMesh().faceLabels();
    const labelList& subAreaFaceLabels = subAreaMesh().faceLabels();

    faceSubToBaseAreaMapPtr_ = new labelList(subAreaFaceLabels.size(), -1);

    labelList& faceSubToBaseAreaMap = *faceSubToBaseAreaMapPtr_;

    Map<label> basePolyFaceHashMap;
    forAll (baseAreaFaceLabels, baseAreaFaceI)
    {
        label basePolyFaceI = baseAreaFaceLabels[baseAreaFaceI];

        if (basePolyFaceI != -1)
        {
            basePolyFaceHashMap.insert(basePolyFaceI, baseAreaFaceI);
        }
    }

    forAll (subAreaFaceLabels, subAreaFaceI)
    {
        label subPolyFaceI = subAreaFaceLabels[subAreaFaceI];
        label basePolyFaceI = faceSubToBaseMap()[subPolyFaceI];

        Map<label>::iterator iter =
            basePolyFaceHashMap.find(basePolyFaceI);

        if (iter != basePolyFaceHashMap.end())
        {
            label baseAreaFaceI = iter();

            faceSubToBaseAreaMap[subAreaFaceI] = baseAreaFaceI;
        }
    }
}

void faSubMesh::makeFaceCurvatures() const
{
    if (debug)
    {
        Info<< "faSubMesh::makeFaceCurvatures() : "
            << "making face curvatures"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (faceCurvaturesPtr_)
    {
        FatalErrorIn("faSubMesh::makeFaceCurvatures()")
            << "face curvatures field already exists"
                << abort(FatalError);
    }

    faceCurvaturesPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "faceCurvatures",
                basePolyMesh().pointsInstance(),
                faMesh::meshSubDir,
                basePolyMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            baseAreaMesh().faceCurvatures()
        );

    areaScalarField& baseK = *faceCurvaturesPtr_;
    scalarField& baseKin = baseK.internalField();
    baseKin = scalarField(baseK.size(), 0);

    const areaScalarField& subK = subAreaMesh().faceCurvatures();
    const scalarField& subKin = subK.internalField();

    scalarField baseKinWeights(baseKin, 0);

    forAll (subKin, subFaceI)
    {
        label baseFaceI = faceSubToBaseAreaMap()[subFaceI];

        baseKin[baseFaceI] += subKin[subFaceI];

// TODO: Use signed distance weights!
        baseKinWeights[baseFaceI] += 1.0;
    }

    baseKin /= baseKinWeights;

    baseK.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& faSubMesh::faceSubToBaseAreaMap() const
{
    if (!faceSubToBaseAreaMapPtr_)
    {
       makeFaceSubToBaseAreaMap();
    }

    return *faceSubToBaseAreaMapPtr_;
}

const areaScalarField& faSubMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        makeFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
