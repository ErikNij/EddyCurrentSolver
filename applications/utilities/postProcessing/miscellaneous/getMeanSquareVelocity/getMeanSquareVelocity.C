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
    getMeanSquareVelocity

Description
    Get mean-square velocty

Author
    Pascal Beckstein

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "OStringStream.H"
#include "demandDrivenData.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    // Add region option
    #include "addRegionOption.H"

    // Add time selector options
    timeSelector::addOptions();

    // Add export file name option
    argList::validOptions.insert("file", "name");

    // Set root case and create time
    #include "setRootCase.H"
    #include "createTime.H"

    // Create time instant list
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

    // Info about what is going on
    Info << "Mean square velocity data" << endl << nl;

    // Open export file if necessary and swap output pointer
    word exportFileName;
    bool exportFileOption = args.optionReadIfPresent<word>("file", exportFileName);

    // Create output file stream and temporary export stream
    OFstream File(exportFileName);
    OStringStream* Export = NULL;

    // Info if export file is beeing used
    if (exportFileOption) { Info << "Export to file: " << exportFileName << endl; }

    // Write info line
    Export = new OStringStream();
    *Export << "# 1:<time>, 2:<Ux^2>, 3:<Uy^2>, 4:<Uz^2>" << endl;
    if (exportFileOption) { File << Export->str().c_str(); } else { Info << "> " << Export->str().c_str(); }
    deleteDemandDrivenData(Export);

    // Loop over all time instances
    forAll(times, timeI)
    {
        // Set time and update mesh
        runTime.setTime(times[timeI], timeI);
        mesh.readUpdate();

        // Read U object
        IOobject UObject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        );

        // Check if U exists
        if (UObject.headerOk())
        {
            mesh.readUpdate();

            // Read velocity field
            volVectorField U(UObject, mesh);

            // Init export variables
            vector  U2 = vector::zero;
            scalar vol = 0.0;

            //      add up the square of the velocity components and then divide by the cell volume
            forAll(mesh.C(), cellI)
            {
                U2[0] += Foam::pow(U[cellI][0],2.0) * mesh.V()[cellI];
                U2[1] += Foam::pow(U[cellI][1],2.0) * mesh.V()[cellI];
                U2[2] += Foam::pow(U[cellI][2],2.0) * mesh.V()[cellI];
                vol   += mesh.V()[cellI];
            }

            // Write data
            Export = new OStringStream();
            *Export << runTime.timeName() << " " << U2[0]/vol << " " << U2[1]/vol << " " << U2[2]/vol << endl;
            if (exportFileOption) { File << Export->str().c_str(); } else { Info << "> " << Export->str().c_str(); }
            deleteDemandDrivenData(Export);

        }
        else
        {
            // Write comment line to show that U is missing
            Export = new OStringStream();
            *Export << "# No velocity data for time " << runTime.timeName() << endl;
            if (exportFileOption) { File << Export->str().c_str(); } else { Info << "> " << Export->str().c_str(); }
            deleteDemandDrivenData(Export);
        }

    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
