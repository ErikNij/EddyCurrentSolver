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

Application
    extractPatchFaceCenters

Description
    Extract patch face centers

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    // Add command line options for time instant selection
    timeSelector::addOptions();

    // Add command line argument to provide a patch name
    argList::validArgs.append("patchName"); // arg #1

    // Add export file name option
    argList::validOptions.insert("file", "name");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList times = timeSelector::select0(runTime, args);

    // Read region name
    Foam::word meshRegionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", meshRegionName);

    // Create Mesh
    fvMesh mesh
    (
        IOobject
        (
            meshRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Loop over all times from instant list
    forAll(times, timeI)
    {

        // Set time and update mesh
        runTime.setTime(times[timeI], timeI);
        mesh.readUpdate();

        // Try to find patch label for given patch name from command line arg #1
        const word exPatchName (IStringStream(args.args()[1])());
        const label exPatchID = mesh.boundaryMesh().findPatchID(exPatchName);

        if (exPatchID > -1)
        {

            // Reference to the extraction data
            const vectorField& exData = mesh.boundaryMesh()[exPatchID].faceCentres();

            // Set file names
            fileName exDir(runTime.timeName()/"rawdata");
            fileName exFile(exDir/"dataCenters."+exPatchName+".dat");

            // Print info
            Info << "Extract face center data to file: " << exFile << endl << endl;

            // Create output directory and stream
            mkDir(exDir);
            OFstream ex(exFile);

            ex << "# Source = Patch, ID: " << exPatchID << ", Name: " << exPatchName << endl;
            ex << "# Time   = " << runTime.timeName() << endl;
            ex << "# Data   = Face center coordinates (1:x 2:y 3:z)" << endl;

            // Extract all data
            forAll(exData, exDataI)
            {

                ex << exData[exDataI].component(0) << " "
                   << exData[exDataI].component(1) << " "
                   << exData[exDataI].component(2) << endl;

            }

        } else {

              FatalError
                  << "A patch called '" << exPatchName << "' does not exist!"
                  << abort(FatalError);

        }
    }

    Info << "End" << endl << endl;

}

// ************************************************************************* //
