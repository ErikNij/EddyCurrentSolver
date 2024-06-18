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

#include "pimpleAppControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::pimpleApp::Control, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool Foam::pimpleApp::Control::criteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict =
        mesh_.solutionDict().solverPerformanceDict();

    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        const label fieldI = applyToField(variableName);
        if (fieldI != -1)
        {
            const List<solverPerformanceData> spd = iter().stream();

            // use either stored residual or read from dict
            scalar residual = GREAT;
            if (residualStorage_[fieldI].stored)
            {
                residual = residualStorage_[fieldI].residual;
                residualStorage_[fieldI].stored = false;
            }
            else
            {
                residual = spd.last().maxInitialResidual();
            }


            checked = true;

            if (storeIni)
            {
                residualControl_[fieldI].initialResidual =
                    spd.first().maxInitialResidual();
            }

            const bool absCheck = residual < residualControl_[fieldI].absTol;
            bool relCheck = false;

            scalar relative = 0.0;
            if (!storeIni)
            {
                const scalar iniRes =
                    residualControl_[fieldI].initialResidual
                  + ROOTVSMALL;

                relative = residual/iniRes;
                relCheck = relative < residualControl_[fieldI].relTol;
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug > 1)
            {
                Info<< algorithmName_ << " loop:" << endl;

                Info<< "    " << variableName
                    << " " << algorithmName_ << " iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldI].initialResidual
                    << ", abs tol = " << residual
                    << " (" << residualControl_[fieldI].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldI].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleApp::Control::Control(fvMesh& mesh, const word& dictName)
:
    pimpleControl(mesh, dictName),
    pRefCell_(0),
    pRefValue_(0.0),
    residualStorage_(residualControl_.size())
{
    // init residual storage
    forAll(residualStorage_, fieldI)
    {
        residualStorage_[fieldI].stored = false;
        residualStorage_[fieldI].residual = GREAT;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleApp::Control::~Control()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pimpleApp::Control::storeResiduals(const word& name)
{
    const dictionary& solverDict =
        mesh_.solutionDict().solverPerformanceDict();

    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();

        // find field name
        if (variableName == name)
        {
            const label fieldI = applyToField(variableName);
            if (fieldI != -1)
            {
                const List<solverPerformanceData> spd = iter().stream();

                // store current residual
                residualStorage_[fieldI].residual =
                    spd.last().maxInitialResidual();

                // indicate that the residual has been stored
                residualStorage_[fieldI].stored = true;

                if (debug > 1)
                {
                    Info<< algorithmName_ << " loop:" << endl;

                    Info<< "    " << variableName
                        << " " << algorithmName_ << " iter " << corr_
                        << ": storing res = "
                        << residualStorage_[fieldI].residual
                        << endl;
                }
            }

            break;
        }
    }
}

void Foam::pimpleApp::Control::skipZeroNonOrtho(const word& name)
{
    // skip if this is the final non-orthogonal iteration
    if (!finalNonOrthogonalIter())
    {
        int nIter = -1;

        const dictionary& solverDict =
            mesh_.solutionDict().solverPerformanceDict();

        forAllConstIter(dictionary, solverDict, iter)
        {
            const word& variableName = iter().keyword();

            // find field name
            if (variableName == name)
            {
                const label fieldI = applyToField(variableName);
                if (fieldI != -1)
                {
                    const List<solverPerformanceData> spd = iter().stream();

                    nIter = spd.last().nIterations();
                }

                break;
            }
        }

        if (nIter == 0)
        {
            if (debug > 1)
            {
                Info<< "    " << name
                    << " " << algorithmName_ << " iter " << corr_
                    << ": nIter = "
                    << nIter
                    << ", skipping further non-orthogonal correctors"
                    << " after the next one"
                    << endl;
            }

            // stop non-orthognal correction after next (final) iteration
            corrNonOrtho_ = nNonOrthCorr_;
        }
    }
}

// ************************************************************************* //
