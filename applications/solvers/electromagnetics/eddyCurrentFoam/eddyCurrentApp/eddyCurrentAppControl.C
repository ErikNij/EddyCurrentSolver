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

#include "eddyCurrentAppControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::eddyCurrentApp::Control, 1);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::eddyCurrentApp::Control::read()
{
    solutionControl::read(false);

    // Read solution controls
    const dictionary& eddyCurrentDict = dict();
    nCorrEDDYCURRENT_ =
        eddyCurrentDict.lookupOrDefault<label>("nCorrectors", 1);
    nSubCorrEDDYCURRENT_ =
        eddyCurrentDict.lookupOrDefault<label>("nSubCorrectors", 2147483647);

    // Manual override for solution direction. For experts only!
    Vector<label>& solutionDir = const_cast<Vector<label>&>(solutionDir_);
    solutionDir =
        eddyCurrentDict.lookupOrDefault<Vector<label> >("solutionDir", solutionDir_);

#ifdef eddyCurrentAppLink_H

    // Read update settings
    emSettings_.enabled =
        emUpdateSettingsDict_.lookupOrDefault<bool>("enabled", true);
    emSettings_.outputTimeIndexCycle =
        emUpdateSettingsDict_.lookupOrDefault<int>("outputTimeIndexCycle", 2147483647);
    emSettings_.timeIndexCycle =
        emUpdateSettingsDict_.lookupOrDefault<int>("timeIndexCycle", 2147483647);
    emSettings_.timeCycle =
        emUpdateSettingsDict_.lookupOrDefault<scalar>("timeCycle", VGREAT);
    emSettings_.relDeltaAmax =
        emUpdateSettingsDict_.lookupOrDefault<scalar>("relDeltaAmax", 0.01);

    if (debug > 1)
    {
        Info<< "eddyCurrentApp::Control::read() : "
            << "emUpdateSettingsDict = " << emUpdateSettingsDict_;
    }

    // Read update data
    emUpdateDataDict_.readIfModified();
    emUpdateData_.update = false;
    emUpdateData_.counter =
        emUpdateDataDict_.lookupOrAddDefault<int>("counter", 0);
    emUpdateData_.outputTimeIndex =
        emUpdateDataDict_.lookupOrAddDefault<int>("outputTimeIndex", 0);
    emUpdateData_.lastTime =
        emUpdateDataDict_.lookupOrAddDefault<scalar>("lastTime", -VGREAT);

    if (debug > 1)
    {
        Info<< "eddyCurrentApp::Control::read() : "
            << "emUpdateDataDict = " << emUpdateDataDict_;
    }


#endif
}


bool Foam::eddyCurrentApp::Control::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    forAll (mesh_, regionI)
    {
        const dictionary& solverDict =
            mesh_[regionI].solutionDict().solverPerformanceDict();

        forAllConstIter(dictionary, solverDict, iter)
        {
            const word& variableName = iter().keyword();
            const label fieldI = applyToField(variableName);
            if (fieldI != -1)
            {
                const List<solverPerformanceData> spd = iter().stream();

                const scalar residual = spd.last().maxInitialResidual();

                checked = true;

                scalar absTol = residualControl_[fieldI].absTol;
                if (tolScales_.found(variableName))
                {
                    absTol = min(tolScales_[variableName] * absTol, 0.1);
                }

                bool absCheck = residual < absTol;
                achieved = achieved && absCheck;

                if (debug > 1)
                {
                    Info<< algorithmName_ << " abs statistics:" << endl;

                    Info<< "    " << variableName
                        << ": res = " << residual
                        << " (" << absTol << ")"
                        << endl;
                }
            }
        }
    }

    return checked && achieved;
};


void Foam::eddyCurrentApp::Control::storePrevIterFields() const
{}


bool Foam::eddyCurrentApp::Control::subCriteriaSatisfied()
{
    // no checks on first sub-iteration - nothing has been calculated yet
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    forAll (mesh_, regionI)
    {
        const dictionary& solverDict =
            mesh_[regionI].solutionDict().solverPerformanceDict();

        forAllConstIter (dictionary, solverDict, iter)
        {
            const word& variableName = iter().keyword();
            const label fieldI = applyToField(variableName);
            if (fieldI != -1)
            {
                const List<solverPerformanceData> spd = iter().stream();

                scalar oldOldResidual = VGREAT;
                if (spd.size() > 2)
                {
                    oldOldResidual = spd[spd.size()-3].maxInitialResidual();
                }

                scalar oldResidual = GREAT;
                if (spd.size() > 1)
                {
                    oldResidual = spd[spd.size()-2].maxInitialResidual();
                }
                const scalar residual = spd.last().maxInitialResidual();

                scalar oldRelative =
                    mag(oldResidual/(oldOldResidual+VSMALL) - 1.0);

                scalar relative =
                    mag(residual/(oldResidual+VSMALL) - 1.0);

                checked = true;

                scalar subRelTol =
                    subScale_ * residualControl_[fieldI].relTol;

                const bool relCheck =
                    (max(oldRelative, relative) <= subRelTol);

                achieved = achieved && relCheck;

                if (debug > 1)
                {
                    Info<< algorithmName_ << " rel statistics:" << endl;

                    Info<< "    " << variableName
                        << ": old rel res = " << oldRelative
                        << ", rel res = " << relative
                        << " (" << subRelTol << "/"
                        << residualControl_[fieldI].relTol << ")"
                        << endl;
                }
            }
        }
    }

    return checked && achieved;
}


#ifdef eddyCurrentAppLink_H

bool Foam::eddyCurrentApp::Control::updateRelDeltaA(label movedRegionI)
{
    emUpdateRelDeltaA_ = false;

    volVectorField& prevC
    (
// TODO: Add member function to objectRegistry to get write access!
        const_cast<volVectorField&>
        (
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("emPrevC")
        )
    );

    if (!emUpdateRelDeltaAInitialized_)
    {
        prevC = mesh_[Region::CONDUCTOR].C();
        prevC.correctBoundaryConditions();

        emUpdateRelDeltaAInitialized_ = true;
    }
    else
    {
        regionVolVectorField C
        (
            IOobject
            (
                "C",
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
                dimLength,
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );

        C[Region::DEFAULT] = mesh_[Region::DEFAULT].C();
        C[movedRegionI] = mesh_[movedRegionI].C();
        C.rmap(movedRegionI);
        C.map(Region::CONDUCTOR);
        C[Region::CONDUCTOR].correctBoundaryConditions();

        scalarField magSqrDeltaC =
            magSqr
            (
                C[Region::CONDUCTOR].internalField()
            - prevC.internalField()
            );

        const volVectorField& jRe =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("jRe");
        const volVectorField& jIm =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("jIm");

        scalarField magSqrj = magSqr(jRe) + magSqr(jIm);

        const volVectorField& ARe =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("ARe");
        const volVectorField& AIm =
            mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("AIm");

        scalarField magSqrA(mesh_[Region::CONDUCTOR].C().size(),0.0);

        if
        (
            mesh_[Region::CONDUCTOR].foundObject<volVectorField>("A0Re")
        && mesh_[Region::CONDUCTOR].foundObject<volVectorField>("A0Im")
        )
        {

            const volVectorField& A0Re =
                mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("A0Re");
            const volVectorField& A0Im =
                mesh_[Region::CONDUCTOR].lookupObject<volVectorField> ("A0Im");

            magSqrA = magSqr(ARe+A0Re) + magSqr(AIm+A0Im);
        }
        else
        {
            magSqrA = magSqr(ARe) + magSqr(AIm);
        }

        volScalarField& relDeltaA
        (
// TODO: Add member function to objectRegistry to get write access!
            const_cast<volScalarField&>
            (
                mesh_[Region::CONDUCTOR].lookupObject<volScalarField> ("emRelDeltaA")
            )
        );

        relDeltaA.internalField() =
            physicalConstant::mu0.value()
            * magSqrDeltaC * sqrt(magSqrj / (magSqrA + SMALL));
        relDeltaA.correctBoundaryConditions();

        if (gMax(relDeltaA) > emSettings_.relDeltaAmax)
        {
            prevC = mesh_[Region::CONDUCTOR].C();
            prevC.correctBoundaryConditions();

            emUpdateRelDeltaA_ = true;
        }
    }

    return emUpdateRelDeltaA_;
}


bool Foam::eddyCurrentApp::Control::update(label movedRegionI)
{
    bool uZeroCounter = updateZeroCounter();
    bool uOutputTimeIndex = updateOutputTimeIndex();
    bool uTimeIndex = updateTimeIndex();
    bool uTime = updateTimeIndex();
    bool uRelDeltaA = false;

    if (movedRegionI != -1)
    {
        // WARNING: Must not be called twice!
        uRelDeltaA = updateRelDeltaA(movedRegionI);
    }

    if (debug > 1)
    {
        Info<< "eddyCurrentApp::Control::update() : "
            << "updateZeroCounter() = " << uZeroCounter << endl;
        Info<< "                                    "
            << "updateOutputTimeIndex() = " << uOutputTimeIndex << endl;
        Info<< "                                    "
            << "updateTimeIndex() = " << uTimeIndex << endl;
        Info<< "                                    "
            << "updateTime() = " << uTime << endl;
        Info<< "                                    "
            << "updateRelDeltaA() = " << uRelDeltaA << endl;
    }

    // Reset outputTimeIndex
    if (uOutputTimeIndex)
    {
        emUpdateData_.outputTimeIndex = 0;
        emUpdateDataDict_.set<int>
        (
            "outputTimeIndex",
            emUpdateData_.outputTimeIndex
        );
    }

    return uZeroCounter
        || uOutputTimeIndex
        || uTimeIndex
        || uTime
        || uRelDeltaA;
}

#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddyCurrentApp::Control::Control
(
    regionFvMesh& mesh,
    const word& dictName
)
:
    solutionControl(mesh[Region::DEFAULT], dictName),
    mesh_(mesh),
    meshIs3D_((mesh.nGeometricD() == 3)),
    solutionDir_(meshIs3D_ ? mesh.geometricD() : -mesh.geometricD()),
    interfacePatchName_
    (
        dict().lookup("interface")
    ),
    interfacePatchLabel_
    (
        mesh_[Region::CONDUCTOR].
            boundaryMesh().findPatchID(interfacePatchName_)
    ),
    nCorrEDDYCURRENT_(-1),
    nSubCorrEDDYCURRENT_(-1),
    switchV_(false),
    tolScales_(),
    subCorr_(0),
    subScale_(1.0)
#ifdef eddyCurrentAppLink_H
    ,
    emUpdateSettingsDict_
    (
        dict().subDict("emUpdate")
    ),
    emUpdateDataDict_
    (
        IOobject
        (
            "emUpdate",
            mesh.time().timeName(),
            "uniform",
            mesh[Region::DEFAULT],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    ),
    emUpdateRelDeltaA_(false),
    emUpdateRelDeltaAInitialized_(false)
#endif
{
    if (interfacePatchLabel_ == -1)
    {
        FatalError
            << "Conductor patch name '"
            << interfaceName()
            << "' does not exist."
            << abort(FatalError);
    }

    read();

    Info<< nl;

    if (residualControl_.empty())
    {
        Info<< algorithmName_ << ": no residual control data found. "
            << "Calculations will employ " << nCorrEDDYCURRENT_
            << " corrector loops"
            << " and " << nSubCorrEDDYCURRENT_
            << " sub-corrector loops" << nl << endl;
    }
    else
    {
        Info<< algorithmName_ << ": max iterations = "
            << nCorrEDDYCURRENT_
            << endl
            << algorithmName_ << ": max sub-iterations = "
            << nSubCorrEDDYCURRENT_
            << endl;
        forAll(residualControl_, i)
        {
            Info<< "    field " << residualControl_[i].name << token::TAB
                << ": rel tol " << residualControl_[i].relTol
                << ", tol " << residualControl_[i].absTol
                << nl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::eddyCurrentApp::Control::loop()
{
    read();

    corr_++;

    if (debug > 1)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    if (corr_ == nCorrEDDYCURRENT_ + 1)
    {
        Info<< algorithmName_ << ": not converged within "
            << nCorrEDDYCURRENT_ << " iterations" << endl;

        corr_ = 0;
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }

    if ((corr_ > 1) && criteriaSatisfied())
    {
        Info<< algorithmName_ << ": converged in " << corr_ - 1
            << " iterations" << endl;

        corr_ = 0;
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }
    else
    {
        Info<< nl;
        Info<< algorithmName_ << ": iteration " << corr_ << endl;
        storePrevIterFields();
    }

    return true;
};


Foam::dictionary Foam::eddyCurrentApp::Control::subDict
(
    label regionI,
    const word& name
) const
{
    dictionary dict = mesh_[regionI].
        solutionDict().subDict("solvers").subDict(name);

    scalar relTol = readScalar(dict.lookup("relTol"));

    scalar relTolScale = VGREAT;
    scalar progressLeft;

    scalar tolerance = readScalar(dict.lookup("tolerance"));
    scalar tolScale = 1.0;
    scalar scaledTol = tolerance;
    if (tolScales_.found(name))
    {
        tolScale = tolScales_[name];
    }

    if ((corr_ == 1) && (subCorr_ == 1))
    {
        progressLeft = 1.0;
    }
    else
    {
        progressLeft = 0.0;

        forAll (mesh_, regioni)
        {
            const dictionary& solverDict =
                mesh_[regionI].solutionDict().solverPerformanceDict();

            forAllConstIter (dictionary, solverDict, iter)
            {
                const word& variableName = iter().keyword();
                const label fieldI = applyToField(variableName);
                if (variableName == name)
                {
                    const List<solverPerformanceData> spd = iter().stream();

                    const scalar residual = spd.last().maxInitialResidual();

                    scalar absTol = tolScale * residualControl_[fieldI].absTol;

                    tolerance = min(tolerance, residualControl_[fieldI].absTol);
                    scaledTol = min(tolScale * tolerance, 0.1);

                    // Residual difference from target residual
                    scalar residualDiff = residual - absTol;

                    // Relative tolerance scaling depending on absolute tol
                    relTolScale = pow(absTol, -1);

                    // Linear convergence progress left
                    progressLeft =
                        max
                        (
                            progressLeft,
                            min(max(residualDiff, 0.0), 1.0)
                        );
                }
            }
        }
    }

    // Relative tolerance decrease function (atan, scaled)
    scalar decRelTol =
        2.0 * atan(subScale_ * relTolScale * progressLeft) / M_PI;

    // Relative tolerance of sub-iteration
    scalar subRelTol = decRelTol * relTol;

    // Overwrite relTol in dictionary
    dict.set<scalar>
    (
        "relTol",
        subRelTol
    );

    // Overwrite tolerance in dictionary
    dict.set<scalar>
    (
        "tolerance",
        scaledTol
    );

    if (debug > 1)
    {
        Info<< algorithmName_ << " sub-loop:" << endl;

        Info<< "    " << name
            << " " << algorithmName_ << " iter dict " << subCorr_
            << ": progress left = " << progressLeft
            << ", rel tol scale = " << relTolScale
            << ", rel tol = " << subRelTol << "/" << relTol
            << ", tol scale = " << tolScale
            << ", tol = " << scaledTol << "/" << tolerance
            << endl;
    }

    return dict;
};


void Foam::eddyCurrentApp::Control::setTolScale
(
    const word& name,
    const scalar& tolScale
)
{
    tolScales_.set(name, tolScale);
};


bool Foam::eddyCurrentApp::Control::subLoop()
{
    subCorr_++;

    if (debug > 1)
    {
        Info<< algorithmName_ << " sub-loop: corr = " << subCorr_ << endl;
    }

    if ((subCorr_ > 1) && criteriaSatisfied())
    {
        subCorr_ = 0;
        subScale_ = 1.0;

        return false;
    }

    if (subCorr_ == nSubCorrEDDYCURRENT_ + 1)
    {
        subCorr_ = 0;
        subScale_ /= 10.0;

        Info<< algorithmName_ << " sub-loop: scale = " << subScale_ << endl;

        return false;
    }

    if ((subCorr_ > 1) && subCriteriaSatisfied() && meshIs3D())
    {
        subCorr_ = 0;
        subScale_ /= 10.0;

        Info<< algorithmName_ << " sub-loop: scale = " << subScale_ << endl;

        return false;
    }
    else
    {
        Info<< algorithmName_ << ": sub-iteration " << subCorr_ << endl;
        storePrevIterFields();
    }

    return true;
};


#ifdef eddyCurrentAppLink_H

bool Foam::eddyCurrentApp::Control::needsUpdate(label movedRegionI)
{
    emUpdateData_.update = false;

    if (emSettings_.enabled)
    {
        // Increase output index if due and write to dictionary
        if (mesh_.time().outputTime())
        {
            emUpdateData_.outputTimeIndex += 1;
            emUpdateDataDict_.set<int>
            (
                "outputTimeIndex",
                emUpdateData_.outputTimeIndex
            );
        }

        emUpdateData_.update = update(movedRegionI);

        if (emUpdateData_.update)
        {
            // Increase update counter and write to dictionary
            emUpdateData_.counter += 1;
            emUpdateDataDict_.set<int>
            (
                "counter",
                emUpdateData_.counter
            );

            // Update last update time and write to dictionary
            emUpdateData_.lastTime = mesh_.time().value();
            emUpdateDataDict_.set<scalar>
            (
                "lastTime",
                emUpdateData_.lastTime
            );
        }
    }

    return emUpdateData_.update;
}

#endif


// ************************************************************************* //
