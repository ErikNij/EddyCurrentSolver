// const dictionary& fs = mesh.solutionDict().subDict("FreeSurface");

// Switch smooth
// (
//     fs.lookup("smooth")
// );

// Switch ddtPhiCorr
// (
//     fs.lookup("ddtPhiCorr")
// );

bool ddtPhiCorr =
    mesh.solutionDict().subDict("PISO")
   .lookupOrDefault<Switch>("ddtPhiCorr", false);

// int nFreeSurfCorr = 1;
// if (fs.found("nFreeSurfCorr"))
// {
//     nFreeSurfCorr = readInt(fs.lookup("nFreeSurfCorr"));
// }
