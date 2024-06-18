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
    getFluidDomainVolume

Description
    Get fluid domain volume

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
    Info << "Fluid domain volume data" << endl << nl;

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
    *Export << "# 1:<time>, 2:<iniFluidVolumePhase1>, 3:<iniFluidVolumePhase2>, 4:<curFluidVolumePhase1>, 5:<curFluidVolumePhase2>" << endl;
    if (exportFileOption) { File << Export->str().c_str(); } else { Info << "> " << Export->str().c_str(); }
    deleteDemandDrivenData(Export);

    // Init/Reset inital fluid volumes
    scalar iniFluidVolumePhase1 = -1.0;
    scalar iniFluidVolumePhase2 = -1.0;
    scalar iniFluidVolumeTotal = -1.0;

    // Loop over all time instances
    forAll(times, timeI)
    {
        // Set time and update mesh
        runTime.setTime(times[timeI], timeI);
        mesh.readUpdate();

        // Read fluid volume dict object
        IOobject fluidVolumeObject
        (
            "fluidVolume",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT
        );

        // Init/Reset current fluid volumes
        scalar curFluidVolumePhase1 = -1.0;
        scalar curFluidVolumePhase2 = -1.0;
        scalar curFluidVolumeTotal = -1.0;

        if (fluidVolumeObject.headerOk())
        {
            // Read fluid volume dictionary
            IOdictionary fluidVolumeDict
            (
                IOobject
                (
                    "fluidVolume",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            // Read initial volume
            dictionary& fluidVolumeInitialSubDict = fluidVolumeDict.subDict("initial");
            iniFluidVolumePhase1 = pTraits<scalar>(fluidVolumeInitialSubDict.lookup("phase1"));
            iniFluidVolumePhase2 = pTraits<scalar>(fluidVolumeInitialSubDict.lookup("phase2"));
            iniFluidVolumeTotal = iniFluidVolumePhase1 + iniFluidVolumePhase2;

            // Read current volume
            dictionary& fluidVolumeCurrentSubDict = fluidVolumeDict.subDict("current");
            curFluidVolumePhase1 = pTraits<scalar>(fluidVolumeCurrentSubDict.lookup("phase1"));
            curFluidVolumePhase2 = pTraits<scalar>(fluidVolumeCurrentSubDict.lookup("phase2"));
            curFluidVolumeTotal = curFluidVolumePhase1 + curFluidVolumePhase2;
        }
        else
        {
            // Read fluid indicator field
            volScalarField alpha1
            (
                IOobject
                (
                    "alpha1",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            // Calculate volume for first time dir and use it as a reference
            if (timeI == 0)
            {
                iniFluidVolumePhase1 = fvc::domainIntegrate(alpha1).value();
                iniFluidVolumePhase2 = fvc::domainIntegrate(1.0-alpha1).value();
                iniFluidVolumeTotal = iniFluidVolumePhase1 + iniFluidVolumePhase2;
            }

            // Calculate current volume
            curFluidVolumePhase1 = fvc::domainIntegrate(alpha1).value();
            curFluidVolumePhase2 = fvc::domainIntegrate(1.0-alpha1).value();
            curFluidVolumeTotal = curFluidVolumePhase1 + curFluidVolumePhase2;

            // Create fluid volume dictionary
            IOdictionary fluidVolumeDict
            (
                IOobject
                (
                    "fluidVolume",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                )
            );

            // Write initial volume
            dictionary fluidVolumeInitialSubDict;
            fluidVolumeInitialSubDict.add<scalar>("phase1", iniFluidVolumePhase1);
            fluidVolumeInitialSubDict.add<scalar>("phase2", iniFluidVolumePhase2);
            fluidVolumeInitialSubDict.add<scalar>("total", iniFluidVolumeTotal);
            fluidVolumeDict.add("initial", fluidVolumeInitialSubDict);

            // Write current volume
            dictionary fluidVolumeCurrentSubDict;
            fluidVolumeCurrentSubDict.add<scalar>("phase1", curFluidVolumePhase1);
            fluidVolumeCurrentSubDict.add<scalar>("phase2", curFluidVolumePhase2);
            fluidVolumeCurrentSubDict.add<scalar>("total", curFluidVolumeTotal);
            fluidVolumeDict.add("current", fluidVolumeCurrentSubDict);

            runTime.writeNow();
        }

        // Write data
        Export = new OStringStream();
        *Export << runTime.timeName() << " "
                << iniFluidVolumePhase1 << " "
                << iniFluidVolumePhase2 << " "
                << curFluidVolumePhase1 << " "
                << curFluidVolumePhase2 << endl;
        if (exportFileOption) { File << Export->str().c_str(); } else { Info << "> " << Export->str().c_str(); }
        deleteDemandDrivenData(Export);
    }

    Info << "\nEnd\n" << endl;

    return(0);
}

// ************************************************************************* //
