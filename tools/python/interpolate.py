#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Extended inverse distance interpolation script
# August 2013
# Dr. Thomas Wondrak (t.wondrak@hzdr.de) [basic program idea/structure]
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

# Additional modules
import time
import argparse
import numpy
import scipy

from scipy.spatial import cKDTree
from math import fmod, ceil

from mpi4py import MPI

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Index partitioning for MPI
def indexPartitioning(imax):

    # Serial mode
    if (mpi_size==0):
        return 0, imax

    # Parallel mode
    else:
        # Local index partition
        di = int(ceil(float(imax) / mpi_size))
        # Calculate range
        i0 = di * mpi_rank
        i1 = min(di*(mpi_rank+1), imax)

        return i0, i1



def getQueryInfo(sourceDataNodesProximityIDs, sourceDataNodesProximityDists):

    # Get maximum number of proximity IDs
    if type(sourceDataNodesProximityIDs) == numpy.ndarray:
        maxValIDs = sourceDataNodesProximityIDs.size
    elif type(sourceDataNodesProximityIDs) == int:
        maxValIDs = 1
    else:
        if (mpi_rank == 0):
            raise TypeError("\nERROR: Variable sourceDataNodesProximityIDs has wrong type!")
        sys.exit(1)

    if type(sourceDataNodesProximityDists) == numpy.ndarray:
        minValDists = min(sourceDataNodesProximityDists)
    elif type(sourceDataNodesProximityDists) == float:
        minValDists = sourceDataNodesProximityDists
    elif type(sourceDataNodesProximityDists) == int:
        minValDists = sourceDataNodesProximityDists
    else:
        if (mpi_rank == 0):
            raise TypeError("\nERROR: Variable sourceDataNodesProximityDists has wrong type!")
        sys.exit(1)

    return minValDists, maxValIDs



# Inverse distance
def inverseDistance(svs, ds, n, c, p, z):

    # Derivative interpolation variables
    tv = scipy.zeros(c)
    ws = 0

    # Loop through proximity
    for i in range(n):
        # Catch n=1
        if n == 1:
            d  = ds
            sv = svs[:]
        else:
            d = ds[i]
            sv = svs[i,:]
        # Replace zero distance occurrences with smallest possible values
        if (d == 0): d = z
        # Calc weights
        w  = d ** -p
        # Add up weights to get weights sum
        ws += w
        # Add up weighted field values
        tv += w * sv

    # Divide field values by total weight
    tv /= ws

    # Return interpolated target field value
    return tv



# Human readable time string
def humanReadableTime(time):

    # Create time string based on size
    if (0 <= time < 1):
        timeString = "< 1s"
    elif (1 <= time < 60):
        timeString = str(int(time)) + "s"
    elif (60 <= time < 3600):
        timeString = str(round(time/60,1)) + "m"
    elif (3600 <= time < 216000):
        timeString = str(round(time/3600,1)) + "h"
    elif (216000 <= time < 1512000):
        timeString = str(round(time/216000,1)) + "d"
    else:
        timeString = str(round(time/1512000,1)) + "w"

    # Retuen time string
    return timeString



# Progressbar
def update_progress(keyname, progress, dt):

    # Set progress bar length
    barLength = 54
    # Convert progress variable to float if necessary
    if (isinstance(progress, int)):
        progress = float(progress)

    # Calculate estimated time from elapsed time and progress
    dte=dt/(progress+1e-20)*(1-progress)

    # Convert times to time strings
    dts = humanReadableTime(dt)
    dtes = humanReadableTime(dte)

    # Modifiy status, border progress and calculate progress bar block count
    if (progress < 0):
        progress = 0
        status = "Halt                                "
    elif (0 <= progress < 1):
        status = "(" + dts + " / est. " + dtes + ")   "
    else:
        progress = 1
        status = "Done in " + dts + "                 \r\n"
    block = int(round(barLength*progress))

    # Set the text string for stdout
    text = "\r- " + keyname + ": [{0}] {1}% {2}".format( "="*block + " "*(barLength-block), round(progress*100,1), status)

    # Write to stdout
    sys.stdout.write(text)
    sys.stdout.flush()


def checkIfEmptyFile(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

# --------------------------------------------------------------------------- #
# --- Argument parsing ------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Description
parser = argparse.ArgumentParser(description="This python script realizes a inverse distance interpolation from data given on nodes (source data) to a different set of arbitrary nodes (target nodes). If necessary the script can be run in parallel with MPI. The result (target data) will include both, target nodes and the interpolated data!")

# Basic optional arguments
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="count")

# Basic interpolation parameters
parser.add_argument("-d", "--dimension", metavar="DIMENSION", type=int, default=3, choices=[2,3], help="Interpolation dimension (default: 3)")
parser.add_argument("-z", "--zero", metavar="ZERO", type=float, default=1e-20, help="Interpolation zero replacement (default: 1e-20)")

# Interpolation optional arguments
parser.add_argument("-Pp", "-p", "--InterpolateExpontent", metavar="P_EXPONENT", type=int, default=3, help="Interpolation: Weighting exponent (default: 3)")
parser.add_argument("-Pn", "-n", "--InterpolatePoints", metavar="P_POINTS", type=int, default=8, help="Interpolation: Number of interpolation points (default: 8)")
parser.add_argument("-Pr", "-r", "--InterpolateRadius", metavar="P_RADIUS", type=float, default=numpy.inf, help="Interpolation: Search radius for interpolation points (default: \"inf\")")

# Smoothing optional arguments
parser.add_argument("-S", "--Smooth", action="store_true", help="Smoothing: Switch (default: false)")
parser.add_argument("-Sp", "--SmoothExpontent", metavar="S_EXPONENT", type=int, default=1, help="Smoothing: Weighting exponent (default: 1)")
parser.add_argument("-Sn", "--SmoothPoints", metavar="S_POINTS", type=int, default=3, help="Smoothing: Number of interpolation points (default: 3)")
parser.add_argument("-Sr", "--SmoothRadius", metavar="S_RADIUS", type=float, default=numpy.inf, help="Smoothing: Search radius for interpolation points (default: \"inf\")")
parser.add_argument("-Sf", "--SmoothFactor", metavar="S_FACTOR", type=float, default=0.5, help="Smoothing: Intensity factor (default: 0.5)")

# File optional parameters
parser.add_argument("-c", "--columns", metavar="COLS", nargs=2, type=int, default=[3,0], help="number of columns containing field data and relative offset (default: 3 0)")
parser.add_argument("-i", "--ignore", metavar="IGNORE", nargs=2, type=str, default=["#","#"], help="ignore line prefix for source and tnodes file (default: # #)")
parser.add_argument("-s", "--skip", metavar="SKIP", nargs=2, type=int, default=[0,0], help="skip line count for source and tnodes file (default: 0 0)")
parser.add_argument("-csd", "--columns-source-data", metavar="COLS-SOURCE-DATA", type=int, default=3, help="Number of columns containing field data (default: 3) [replaces -c]")
parser.add_argument("-cso", "--columns-source-offset", metavar="COLS-SOURCE-OFFSET", type=int, help="Absolute columns offset for field data counted from (default: DIM) [replaces -c]")
parser.add_argument("-cno", "--columns-tnodes-offset", metavar="COLS-TNODES-OFFSET", type=int, default=0, help="Absolute columns offset for node data counted from (default: 0)")

# Positional arguments
parser.add_argument("source", metavar="SOURCE-FILE", nargs="?", default="source.dat", help="source data file (default: \"source.dat\")")
parser.add_argument("tnodes", metavar="TNODES-FILE", nargs="?", default="tnodes.dat", help="target nodes file (default: \"tnodes.dat\")")
parser.add_argument("target", metavar="TARGET-FILE", nargs="?", default="target.dat", help="target data file (default: \"target.dat\")")

# Commence parsing
cargs = parser.parse_args()

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Basic parameters ++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

verb = cargs.verbose



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Inverse distance parameters +++++++++++++++++++++++++++++++++++++++++++ #

# General parameters
d = cargs.dimension                # Dimension
z = cargs.zero                     # Zero replacement

# Interpolation parameters
P = True                           # Interpolation switch
Pp = cargs.InterpolateExpontent    # Interpolation exponent
Pn = cargs.InterpolatePoints       # Interpolation points
Pr = cargs.InterpolateRadius       # Interpolation radius

# Smoothing parameters
S = cargs.Smooth                   # Smoothing switch
if not S: S = False
Sp = cargs.SmoothExpontent         # Smoothing exponent
Sn = cargs.SmoothPoints            # Smoothing points
Sr = cargs.SmoothRadius            # Smoothing radius
Sf = cargs.SmoothFactor            # Smoothing factor



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ File and column parameters ++++++++++++++++++++++++++++++++++++++++++++++ #

c = cargs.columns[0]            # Inverse distance field data columns count
o = cargs.columns[1]            # Inverse distance field data columns offset

# Replacement options (this is only to keep the -c option for backward compatibility)
csd = c
cso = o + d # Default for cso
cno = cargs.columns_tnodes_offset
if not cargs.columns_source_data == None: csd = cargs.columns_source_data
if not cargs.columns_source_offset == None: cso = cargs.columns_source_offset

# TODO: Implement usage of cno


# Source parameters
sourceFile = cargs.source          # Source file name
sourceFileSkip = cargs.skip[0]     # Source file header line skip value
sourceFileIgnore = cargs.ignore[0] # Source file header/comment line prefix

# Target parameters
tnodesFile = cargs.tnodes          # Target nodes file name
tnodesFileSkip = cargs.skip[1]     # Target nodes file header line skip value
tnodesFileIgnore = cargs.ignore[1] # Target nodes file header/comment line prefix
targetFile = cargs.target          # Target file name



# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ MPI init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Initialize MPI and say hello to the world communicator
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Prepare file and column data ++++++++++++++++++++++++++++++++++++++++++ #

# Set columns for nodes and field values by dimension
cols = range(cso+csd)
colsNodes = cols[:d]
colsFields = cols[cso:cso+csd]

# Read source data, create nodes and field parts, count rows
if (checkIfEmptyFile(sourceFile)):
    sourceData  = scipy.loadtxt(sourceFile, skiprows=sourceFileSkip, comments=sourceFileIgnore)
else:
    if (mpi_rank == 0): raise IOError("ERROR: Empty source file!")
sourceDataNodes = sourceData[:, colsNodes]
sourceDataNodesCount = sourceDataNodes.shape[0]
sourceDataFieldValues = sourceData[:, colsFields]

# Read target nodes, count them and init target field part
if (checkIfEmptyFile(tnodesFile)):
    targetDataNodes = scipy.loadtxt(tnodesFile, skiprows=tnodesFileSkip, comments=tnodesFileIgnore)
else:
    if (mpi_rank == 0): raise IOError("ERROR: Empty target node file!")
targetDataNodesCount = targetDataNodes.shape[0]
targetDataFieldValues = scipy.zeros((targetDataNodesCount, sourceDataFieldValues.shape[1]))



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Settings overview +++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# TODO: Implement info overview



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Basic checks ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Check source node count
if (sourceDataNodes.shape[0] == 0):
    if (mpi_rank == 0): raise ValueError("ERROR: Source node count is zero!")

# Check target node count
if (targetDataNodes.shape[0] == 0):
    if (mpi_rank == 0): raise ValueError("ERROR: Target node count is zero!")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Inverse distance interpolation ++++++++++++++++++++++++++++++++++++++++ #

if P:

    # Built up the search tree
    if (verb > 1 and mpi_rank == 0):
        print("Building up k-d-tree for interpolation...")
    sourceDataNodesTree = cKDTree(sourceDataNodes)

    # Get local target node range from MPI index partitioning
    [i0,i1] = indexPartitioning(targetDataNodesCount)

    # Get ready and init average proximity distribution matrix on verbose
    if (mpi_rank == 0 and verb > 1):
        print("Interpolating now...")

    # Init progress bar
    timeStart = time.time()
    if (mpi_rank == 0 and verb > 0):
        update_progress("Interpolate",0,0)

    # Loop over all local target nodes
    for i in range(i0,i1):

        # Update time
        timeElapsed = time.time() - timeStart

        # Progress bar update
        if (mpi_rank == 0 and verb > 0):
            i_progress = (i+1-i0)/float((i1-i0))
            if (int(fmod(i_progress*1000,1)) == 0):
                update_progress("Interpolate",i_progress,timeElapsed)

        # Search source node neighbours for current target node and store their euclidian distance and IDs
        sourceDataNodesProximityDists, sourceDataNodesProximityIDs = sourceDataNodesTree.query(targetDataNodes[i], k=Pn, p=2, distance_upper_bound=Pr)

        # Check proximity IDs and stop if not enough points were found
        sourceMinValDists, sourceMaxValIDs = getQueryInfo(sourceDataNodesProximityIDs, sourceDataNodesProximityDists)
        if not (sourceMaxValIDs > 0 or sourceMinValDists < numpy.inf):
            if (mpi_rank == 0):
                if (verb > 0): update_progress("Interpolate",-1,timeElapsed)
                raise ValueError("\nERROR: Tree query failure! Try to increase search radius or decrease \nnumber of interpolation points!")

        # Inverse distance calculation
        targetDataFieldValues[i,:] = inverseDistance(sourceDataFieldValues[sourceDataNodesProximityIDs,:], sourceDataNodesProximityDists, sourceMaxValIDs, csd, Pp, z)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Inverse distance smoothing ++++++++++++++++++++++++++++++++++++++++++++ #

if S:

    # Built up the search tree for smoothing
    if (verb > 1 and mpi_rank == 0):
        print("Building up k-d-tree for smoothing...")
    smoothDataNodesTree = cKDTree(targetDataNodes)

    # Get local target node range from MPI index partitioning
    [i0,i1] = indexPartitioning(targetDataNodesCount)

    # Get ready and init average proximity distribution matrix on verbose
    if (mpi_rank == 0 and verb > 1):
        print("Smoothing now...")

    # Init progress bar
    timeStart = time.time()
    if (mpi_rank == 0 and verb > 0):
        update_progress("Smooth",0,0)

    # Loop over all local target nodes
    for i in range(i0,i1):

        # Update time
        timeElapsed = time.time() - timeStart

        # Progress bar update
        if (mpi_rank == 0 and verb > 0):
            i_progress = (i+1-i0)/float((i1-i0))
            if (int(fmod(i_progress*1000,1)) == 0):
                update_progress("Smooth",i_progress,timeElapsed)

        # Search source node neighbours for current target node and store their euclidian distance and IDs
        smoothDataNodesProximityDists, smoothDataNodesProximityIDs = smoothDataNodesTree.query(targetDataNodes[i], k=Sn+1, p=2, distance_upper_bound=Sr)

        # Remove owner
        if (smoothDataNodesProximityDists.size == 2):
            smoothDataNodesProximityDists = float(smoothDataNodesProximityDists[1])
        else:
            smoothDataNodesProximityDists = smoothDataNodesProximityDists[1:smoothDataNodesProximityDists.size]

        if (smoothDataNodesProximityIDs.size == 2):
            smoothDataNodesProximityIDs = int(smoothDataNodesProximityIDs[1])
        else:
            smoothDataNodesProximityIDs = smoothDataNodesProximityIDs[1:smoothDataNodesProximityIDs.size]

        # Check proximity IDs and stop if not enough points were found
        smoothMinValDists, smoothMaxValIDs = getQueryInfo(smoothDataNodesProximityIDs, smoothDataNodesProximityDists)
        if not (smoothMaxValIDs > 0 or smoothMinValDists < numpy.inf):
            if (mpi_rank == 0):
                if (verb > 0): update_progress("Smooth",-1,timeElapsed)
                raise ValueError("\nERROR: Tree query failure! Try to increase search radius or decrease \nnumber of interpolation points!")

        # Inverse distance calculation
        targetDataFieldValues[i,:] *= (1.-Sf)
        targetDataFieldValues[i,:] += Sf * inverseDistance(targetDataFieldValues[smoothDataNodesProximityIDs,:], smoothDataNodesProximityDists, smoothMaxValIDs, csd, Sp, z)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Local output data  ++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Concatenate local target data
targetDataPart = scipy.concatenate((targetDataNodes, targetDataFieldValues),1)

# Save/Export local target data file
scipy.savetxt(targetFile + str(mpi_rank), targetDataPart[i0:i1,:])



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ MPI barrier  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Barrier MPI
mpi_comm.Barrier()



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Global output data  +++++++++++++++++++++++++++++++++++++++++++++++++++ #

# Combine local target data files to one global file
if (mpi_rank == 0 and mpi_size != 0):

    # Generate list with all local target file parts
    targetFileParts=list(targetFile + str(rank) for rank in range(mpi_size))

    # Read and concatenate local parts
    targetData = "".join([open(targetFilePart).read() for targetFilePart in targetFileParts])

    # Save/Export global target data file
    open(targetFile, "w").write(targetData)

    # Remove local target data files
    for targetFilePart in targetFileParts:
        os.remove(targetFilePart)

# --------------------------------------------------------------------------- #
# --- Finish and exit ------------------------------------------------------- #
# --------------------------------------------------------------------------- #

MPI.Finalize()
sys.exit(0)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
