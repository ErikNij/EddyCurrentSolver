#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Field update script
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# TODO [High]: Extend/Modify interpolate.py to allow list of tnodes/target pairs
#              to be processed in a row, while cKDtree is only created once

# TODO [Medium]: Add foam application/utility to extract patch face centers
#                for all data-patches at once (e.g. extractAllPatchFaceCentres)

# TODO [Low]: Allow seperate files as interpolatation source !!!

# TODO [Low]: Interpolation directly in python (no subProc)

# TODO [Low]: Fix output for subprocesses (no verbose should silence them)

# TODO [Low]: Check if all tools are available and exit if not

# TODO [Low]: Add program description

# TODO [Low]: Rework function descriptions

# TODO [Low]: Reorganize some functions in modules

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

# pyFoam modules
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

# foamTools modules
from foamTools import subProc
from foamTools.fileEdit import replacePatternByMatch
import foamTools.classInfo as foamClassInfo
import foamTools.dictHelper as foamDict
from foamTools.fieldWriteFromRaw import fieldWriteFromRaw

# Additional modules
import glob
import shutil
import fileinput
import re
import argparse
import subprocess
import numpy

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def verbStr(vl):

    # return : verb string (e.g. -vvv)
    #
    # vl     : verb level of type int

    # TODO [Low]: Improve descriptions

    if not vl == None:
        if (vl > 0):
            return "-" + "v"*vl + " "

        if (vl < 0):
            raise ValueError("verbose level must not be negative")
        else:
            return " "
    else:
        return " "



def foamTimeSubDirFilePath(sub, pre, t="", p="", f=""):

    # return : path to raw file
    #
    # sub    : sub directory name of time folder
    # pre    : data file prefix
    # t      : data file type (internalField, boundary)
    # p      : patch name if t=boundary
    # f      : field name

    # TODO [Low]: Improve descriptions

    path = foamTimePath + "/" + sub + "/" + pre

    if sub==foamTimeRawSubDir:
        if not t=="internalField":
            if t=="boundary":
                path += "." + p

            else:
                raise ValueError("the type " + t + "has to be either \"interalField\" or \"boundary\"")

    if not f=="":
        path += "." + f

    if any(sub in s for s in [foamTimeRawSubDir, foamTimeIntSubDir]):
        path += ".dat"

    elif sub==foamTimeLogSubDir:
        path += ".log"

    else:
        raise ValueError("the time sub dir " + sub + "has to be one of " + str([foamTimeRawSubDir, foamTimeIntSubDir, foamTimeLogSubDir]))

    return path

def  foamTimeRawFilePath(pre, t="", p="", f=""):
    return foamTimeSubDirFilePath(foamTimeRawSubDir, pre, t, p, f)

def  foamTimeIntFilePath(pre, t="", p="", f=""):
    return foamTimeSubDirFilePath(foamTimeIntSubDir, pre, t, p, f)

def  foamTimeLogFilePath(pre, t="", p="", f=""):
    return foamTimeSubDirFilePath(foamTimeLogSubDir, pre, t, p, f)



def foamGetPatchListFromBoudaryDict(boundaryDict):

    """
    Read all patches and their base types directly from the boundary
    dictionary of the current polyMesh.
    @return: (patch list, "patch name <-> type mapping" as list in
             same order as the patch list)
    @param boundaryDict: Currently valid boundary dictionary
    """

    # TODO [Low]: Intercept error exceptions with useful descriptions

    if not type(boundaryDict)== ParsedParameterFile:
        raise TypeError(str(boundaryDict) + " is no valid dictionary (ParsedParameterFile)")

    # Read and store patch types
    pl = []
    ptl = []
    for p in range(0,len(boundaryDict),2):
        n = boundaryDict[p]
        t = foamDict.readStr(boundaryDict[p+1]["type"])
        pl.append(boundaryDict[p])
        ptl.append((n, t))

    return pl, ptl



def foamGetUpdateMethodFromFieldUpdateDict(dictionary, argMethod=None):

    """
    Function to determine the update method.
    Depending on the exsitance of argMethod, the method from the fieldUpdateDict
    may be overwritten.
    @return: update method
    @param dictionary : field update dictionary
    @param argMethod: method name read from command line
    """

    if not type(dictionary)== ParsedParameterFile:
        raise TypeError(str(dictionary) + " is no valid dictionary (ParsedParameterFile)")

    # update method from dictionary if no optional argument was given
    if argMethod==None:
        um = foamDict.readStr(dictionary["updateMethod"])

    # Update method from command line option if present
    else:
        um = argMethod

    return um



def foamGetUpdateMethodDictFromFieldUpdateDict(dictionary, updateMethod):

    """
    Function to determine the update method"s dictionary.
    @return: update method"s dictionary
    @param dictionary : field update dictionary
    @param updateMethod: method name
    """

    if foamDict.tryKey(dictionary, updateMethod):
        return foamDict.readSubDict(dictionary[updateMethod])
    else:
        raise KeyError("the dictionary for selected update method "" + updateMethod + "" is missing")



def foamGetFieldListFromFieldUpdateDict(dictionary):

    """
    Function to return a list of all fields and their own sub-dictionaries
    from the main field update dictionary.
    @return: (field list, field dictionary list in same order)
    @param dictionary: field update dictionary
    """

    if not type(dictionary)== ParsedParameterFile:
        raise TypeError(str(dictionary) + " is no valid dictionary (ParsedParameterFile)")

    if not foamDict.tryKey(dictionary, "fields"):
        raise KeyError("the \"fields\" dictionary is missing")

    fl = []
    fdl = []
    for f in foamDict.readSubDict(dictionary["fields"]):
        fl.append(f)
        fdl.append(foamDict.readSubDict(dictionary["fields"][f]))

    return fl, fdl



def foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict(key, sub=""):

    """
    Small helper function to get key values from the parent field as default
    if specific keys or sub-dictionaries are not given in each single patch.
    @param key: Key to search in patch dict of fieldUpdateDict
    @param sub: Sub-dictionary of patch/field to search the key in
    """

    # TODO [Low]: Intercept error exceptions with useful descriptions

    if foamDict.tryKey(fieldDictList[f], "patches"):
        if foamDict.tryKey(fieldDictList[f]["patches"], patch):
            if sub=="":
                if foamDict.tryKey(fieldDictList[f]["patches"][patch], key):
                    return fieldDictList[f]["patches"][patch][key]
            else:
                if foamDict.tryKey(fieldDictList[f]["patches"][patch], sub):
                    if foamDict.tryKey(fieldDictList[f]["patches"][patch][sub],key):
                        return fieldDictList[f]["patches"][patch][sub][key]

    else:
        if sub=="":
            return fieldDictList[f][key]
        else:
            return fieldDictList[f][sub][key]

# --------------------------------------------------------------------------- #
# --- Argument parsing ------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Description
parser = argparse.ArgumentParser(description="")

# Optional arguments
parser.add_argument("-c", "--clear",
                    help="(create and) clear subdirectories (rawdata, intdata, log) of the time folder given with TIME before applying the field update method", action="store_true")
parser.add_argument("-C", "--clean",
                    help="(create and) clear subdirectories (rawdata, intdata, log) of the time folder and exit", action="store_true")
parser.add_argument("-m", "--method", metavar="METHOD",
                    help="overwrite name of update method in", action="store",
                    type=str)
parser.add_argument("-v", "--verbose",
                    help="increase output verbosity", action="count")

# Positional arguments
parser.add_argument("time", metavar="TIME",
                    help="name of time directory to update (default: 0)", action="store",
                    type=str, default="0")

# Commence parsing
cargs = parser.parse_args()

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Argument parsing
clear = cargs.clear
clean = cargs.clean
if clean: clear = True
method = cargs.method
verb = cargs.verbose
if verb == None: verb = 0
time = cargs.time

# methods
methodValidResultTypes = ["foamdata", "rawdata", "intdata"]

# foam
foamCasePath = os.getcwd()

foamTimePath = time
foamTimeFaMeshPath = foamTimePath + "/faMesh"
foamTimePolyMeshPath = foamTimePath + "/polyMesh"
foamTimeBoundaryDictPath = foamTimePolyMeshPath + "/boundary"

foamTimeRawSubDir = "rawdata"
foamTimeRawPath = foamTimePath + "/" + foamTimeRawSubDir
foamTimeIntSubDir = "intdata"
foamTimeIntPath = foamTimePath + "/" + foamTimeIntSubDir
foamTimeLogSubDir = "log"
foamTimeLogPath = foamTimePath + "/" + foamTimeLogSubDir
foamTimeSubPathList = [foamTimeRawPath, foamTimeIntPath, foamTimeLogPath]

foamConstantPath = "constant"
foamConstantFaMeshPath = foamConstantPath + "/faMesh"
foamConstantPolyMeshPath = foamConstantPath + "/polyMesh"
foamConstantBoundaryDictPath = foamConstantPolyMeshPath + "/boundary"

foamSystemPath = "system"

foamToolsPath = "tools"
foamToolsFieldUpdateDictPath = foamToolsPath + "/fieldUpdateDict"

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ print(basic info ++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

print()
if (verb > 0): print("time = " + time)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ parse boundary dictionary to get all patch information ++++++++++++++++ #

# Find current boundary file
if os.path.exists(foamTimeBoundaryDictPath):
    boundaryDictPath = foamConstantBoundaryDictPath
else:
    boundaryDictPath = foamConstantBoundaryDictPath

# Parse boundary dictionary
if (verb > 0): print("read boundary dictionary")
boundaryDict = ParsedParameterFile(boundaryDictPath, boundaryDict=True)

# Extract patch-types mapping as list
if (verb > 0): print("create patch list")
patchList, patchTypeList = foamGetPatchListFromBoudaryDict(boundaryDict)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ parse field update dictionary +++++++++++++++++++++++++++++++++++++++++ #

# Parse field update dictionary
if (verb > 0): print("read field update dictionary")
fieldUpdateDict = ParsedParameterFile(foamToolsFieldUpdateDictPath)
# FIXME: Why is macro expansion not working correctly?
fieldUpdateDict = ParsedParameterFile(foamToolsFieldUpdateDictPath, doMacroExpansion=True)

# Extract field update method
if (verb > 0): print("set update method")
updateMethod = foamGetUpdateMethodFromFieldUpdateDict(fieldUpdateDict, method)

# Extract all field names and parse their sub-dictionaries
if (verb > 0): print("create field list and extract their dictionaries")
fieldList, fieldDictList = foamGetFieldListFromFieldUpdateDict(fieldUpdateDict)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Clean and prepare current time directory ++++++++++++++++++++++++++++++ #

# Check if foam time directory exists
if not (os.path.isdir(time)):
    raise IOError("time directory "" + time + "" does not exist")

# Create foam time directory structure
if (verb > 1): print("create time directory structure")
for dirPath in foamTimeSubPathList:
    if not (os.path.isdir(dirPath)):
        if (verb > 1): print("create " + dirPath)
        os.mkdir(dirPath)

if clear:
    # Raw/Int data and log files
    if (verb > 0): print("delete old files in time sub-directories " + str(foamTimeSubPathList))
    for dirPath in foamTimeSubPathList:
        for path in glob.glob(dirPath + "/" + "*"):
            if not os.path.isdir(path):
                if os.path.exists(path): os.remove(path)

# Stop if we would like to clean only
if clean:
    print()
    sys.exit(0)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Pre-processing depending on selected update method ++++++++++++++++++++ #

# Process dummy process aliases
if updateMethod=="foam":

    updateMethodInternal = "dummy"
    methodResultType = "foamdata"

elif updateMethod=="raw":

    updateMethodInternal = "dummy"
    methodResultType = "rawdata"

elif updateMethod=="int":

    updateMethodInternal = "dummy"
    methodResultType = "intdata"

else:
    updateMethodInternal = updateMethod

    # Print info
    if (verb > 0): print("select update method " +  updateMethod)

    # Extract field update method dictionary
    if (verb > 0): print("get update method dictionary")
    updateMethodDict = foamGetUpdateMethodDictFromFieldUpdateDict(fieldUpdateDict, updateMethod)

    # Parse method"s settings dictionary
    updateMethodSettingsDict = foamDict.readSubDict(foamDict.tryMandatory(updateMethodDict, "settings"))

    # Read method"s result type
    methodResultType = foamDict.readStr(foamDict.tryMandatory(updateMethodDict, "type"))
    if not any(methodResultType in t for t in methodValidResultTypes):
        raise ValueError("result type of selected update method is \"" + methodResultType + "\", but must be one of \"" + str(methodValidResultTypes))

# Set interpolation switch to on if we are dealing with data which needs interpoltion
methodInterpolate = False
if methodResultType=="intdata": methodInterpolate = True

methodCreateFromRaw = False
if methodResultType=="rawdata" or methodResultType=="intdata": methodCreateFromRaw = True



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Extract points/centers from foam if an interpolation is pending +++++++ #

if methodInterpolate:

    # Extract cell centers of volume mesh
    if (verb > 0): print("extract cell centers")
    exe = "extractCellCenters"
    arg = "-time " + time
    log = foamTimeLogFilePath("extractCellCenters")
    subProc.run(exe, arg, log, verb)

    ## Loop over all data patches
    for p, patch in enumerate(patchList):

        # Read field information
        patchBaseType = patchTypeList[p][1]
        patchPrimitiveType, patchValued = foamClassInfo.patchBaseTypeToPrimitiveValued(patchBaseType)

        # Interpolate only data-valued patches
        if patchValued:

            # Extract patch face center data
            if (verb > 0): print("extract face centers of patch " + patch)
            exe = "extractPatchFaceCenters"
            arg = patch + " -time " + time
            log = foamTimeLogFilePath("extractPatchFaceCenters", p=patch)
            subProc.run(exe, arg, log, verb)

            # Extract patch point data
            if (verb > 0): print("extract points of patch " + patch)
            exe = "extractPatchPoints"
            arg = patch + " -time " + time
            log = foamTimeLogFilePath("extractPatchPoints", p=patch)
            subProc.run(exe, arg, log, verb)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Process update method +++++++++++++++++++++++++++++++++++++++++++++++++ #

if (updateMethodInternal == "dummy"):
    pass

elif (updateMethodInternal == "copyFoam"):

    # Read basic settings from update method dictionary
    copyFoamDataPath = foamDict.readStr(foamDict.tryMandatory(updateMethodSettingsDict, "path"))

    # Copy foam file
    for field in fieldList:
        shutil.copy(copyFoamDataPath + "/" + field, foamTimePath)

elif (updateMethodInternal == "copyRaw"):

    # Read basic settings from update method dictionary
    copyRawDataPath = foamDict.readStr(foamDict.tryMandatory(updateMethodSettingsDict, "path"))

    # Copy all data of given path
    for field in fieldList:
        for path in glob.glob(copyRawDataPath + "/" + field + "*"):
            shutil.copy(path, foamTimeRawPath)

elif (updateMethodInternal == "comsolLorentzForceSurface"):

    # Name of comsol binary
    comsolBin = "comsol"

    # Read basic settings from update method dictionary
    comsolWorkPath = foamDict.readStr(foamDict.tryMandatory(updateMethodSettingsDict, "root"))
    comsolRecompile = foamDict.readBool(foamDict.tryOrDefault(updateMethodSettingsDict, "recompile",True))

    # Settings for temporary files and logging
    comsolTmpPath = comsolWorkPath + "/tmp"
    comsolLogPath = comsolWorkPath + "/log"
    comsolSubPathList = [comsolTmpPath, comsolLogPath]

    comsolSurfaceDatFile = "surface.dat"
    comsolSurfaceDatPath = comsolTmpPath + "/" + comsolSurfaceDatFile
    comsolExportDatFile = "export.dat"
    comsolExportDatPath = comsolTmpPath + "/" + comsolExportDatFile

    # Read model name to derive file names
    comsolModelName = foamDict.readStr(foamDict.tryMandatory(updateMethodSettingsDict, "name"))
    comsolModelName = re.sub("\.[^.]*$","",comsolModelName)

    # Read java file path
    comsolJavaFile = comsolModelName + ".java"
    comsolJavaPath = comsolWorkPath + "/" + comsolJavaFile
    comsolJavaPathFull = foamCasePath + "/" + comsolJavaPath

    # Relate class file to java file
    comsolClassFile = comsolModelName + ".class"
    comsolClassPath = comsolWorkPath + "/" + comsolClassFile
    comsolClassPathFull = foamCasePath + "/" + comsolClassPath

    # Parse surface dictionary
    comsolSurfaceDict = foamDict.readSubDict(updateMethodSettingsDict["surface"])

    comsolSurfacePatch = foamDict.readStr(foamDict.tryMandatory(comsolSurfaceDict, "patch"))
    comsolSurfaceSource = foamDict.readStr(foamDict.tryMandatory(comsolSurfaceDict, "source"))
    comsolSurfaceSourcePath = foamTimeRawFilePath(comsolSurfaceSource, t="boundary", p=comsolSurfacePatch)
    comsolSurfaceComments = foamDict.readStr(foamDict.tryMandatory(comsolSurfaceDict, "comments"))

    # Copy result settings from method
    comsolResultType = methodResultType
    comsolResult = foamDict.readStr(foamDict.tryMandatory(updateMethodSettingsDict, "result"))
    comsolResultPath = foamTimeSubDirFilePath(comsolResultType, comsolResult)


    # --- prepare comsol directory ------------------------------------------ #

    # Print info
    if (verb > 0): print("prepare comsol directory structure")

    # Create comsol time directory structure if necessary
    for dirPath in comsolSubPathList:
        if not os.path.isdir(dirPath):
            if (verb > 0): print("create " + dirPath + " in " + time)
            os.mkdir(dirPath)


    # --- prepare surface data ---------------------------------------------- #

    # Print info
    if (verb > 0): print("prepare surface data")

    # Check if surface file is present
    if not os.path.exists(comsolSurfaceSourcePath):
        raise IOError("surface file \"" + comsolSurfaceSourcePath + "\" does not exist")

    # Copy raw surface file for comsol
    if (verb > 0): print("copy and prepare point data file of free surface")
    shutil.copy(comsolSurfaceSourcePath, comsolSurfaceDatPath)
    if not os.path.exists(comsolSurfaceDatPath):
        raise IOError("comsol surface data file \"" + comsolSurfaceDatPath + "\" is no file")

    # Replace comments prefix in surface file
    replacePatternByMatch(comsolSurfaceDatPath, comsolSurfaceComments,"%")


    # --- comsol simulation and field data export --------------------------- #

    # Check if comsol java setup file exists and exit if not
    if not os.path.exists(comsolJavaPath):
        raise IOError("comsol java file \"" + comsolJavaPath + "\" is missing")

    # Recompile comsol java setup file to get a class file if necessary
    if comsolRecompile or not os.path.exists(comsolClassPath):
        if (verb > 0): print("recompile comsol setup")
        exe = comsolBin
        arg = "compile" + " " + comsolJavaPathFull
        log = comsolLogPath + "/compile"
        subProc.run(exe,arg,log,verb)

    # Start comsol
    # HINT: Comsol needs to get started in the directory
    #       where the class file is stored: so we have to
    #       change our working directory, temporarily
    if (verb > 0): print("start comsol simulation")
    os.chdir(comsolWorkPath)
    if True:
        # Run simulation of the setup
        exe = comsolBin
        arg = "batch -inputfile" + " " + comsolClassFile
        log = comsolLogPath + "/batch"
        subProc.run(exe,arg,vrb=verb)
    os.chdir(foamCasePath)

    # Check if comsol export data file exists and exit if not
    if not os.path.exists(comsolExportDatPath):
        raise IOError("comsol export data file \"" + comsolExportDatPath + "\" is missing")


    # --- process results --------------------------------------------------- #

    # Copy epxported data file depending on its type raw/int to time subdir
    if (verb > 0): print("copy comsol results")
    shutil.copy(comsolExportDatPath, comsolResultPath)



else: raise ValueError("update method " + updateMethod + " is not supported")




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# +++ Post-processing (interpolation, field creation) +++++++++++++++++++++++ #

# Loop over all given fields
for f, field in enumerate(fieldList):

    # Read field information
    fieldClass = foamDict.readStr(foamDict.tryMandatory(fieldDictList[f], "class"))
    fieldAtomicType, fieldAtomicSize = foamClassInfo.fieldClassToAtomicTypeSize(fieldClass)
    fieldDimension = foamDict.readDim(foamDict.tryMandatory(fieldDictList[f], "dimension"))

    # Read interpolation switch from method"s settings
    fieldInterpolate = methodInterpolate


    # --- interpolate raw data ------------------------------------------ #

    # Interpolate if necessary
    if fieldInterpolate:

        # Print info
        if (verb > 0): print("interpolate data for field " + field + " on cell centers")

        # Parse interpolation data dictionary for internal field
        intDataDict = foamDict.readSubDict(foamDict.tryMandatory(fieldDictList[f], "intdata"))

        ## Read interpolation source information
        source = foamDict.readStr(foamDict.tryMandatory(intDataDict, "source"))
        sourceColumns = fieldAtomicSize
        sourceOffset = foamDict.readInt(foamDict.tryMandatory(intDataDict, "sourceOffset"))
        sourceComments = foamDict.readStr(foamDict.tryMandatory(intDataDict, "sourceComments"))
        sourcePath = foamTimeIntPath + "/" + source + ".dat"

        # Check interpolation source data of internalField
        if not os.path.exists(sourcePath):
            raise IOError("interpolation source data file \"" + sourcePath + "\" for field " + field + " is missing")

        # Read interpolation target node information
        tnodes = foamDict.readStr(foamDict.tryMandatory(intDataDict, "tnodes"))
        tnodesOffset = foamDict.readInt(foamDict.tryMandatory(intDataDict, "tnodesOffset"))
        tnodesComments = foamDict.readStr(foamDict.tryMandatory(intDataDict, "tnodesComments"))
        tnodesPath = foamTimeRawPath + "/" + tnodes + ".internal" + ".dat"

        # Check interpolation target node data of internalField
        if not os.path.exists(tnodesPath):
            raise IOError("interpolation target node data file \"" + tnodesPath + "\" for field " + field + " is missing")

        # Assemble target information
        target = foamDict.readStr(foamDict.tryOrDefault(intDataDict, "target",field))
        targetPath = foamTimeRawPath + "/" + target + "." + "internal" + ".dat"

        # Parse interpolation data settings sub-dictionary for internal field
        intDataSettingsDict = foamDict.readSubDict(foamDict.tryMandatory(intDataDict, "settings"))

        # Read settings
        settingsInterpolateExponent = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Pp",3))
        settingsInterpolatePoints = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Pn",8))
        settingsInterpolateRadius = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Pr",numpy.inf))
        settingsSmooth = foamDict.readBool(foamDict.tryOrDefault(intDataSettingsDict, "S",False))
        settingsSmoothStr = ""
        if settingsSmooth: settingsSmoothStr = "-S"
        settingsSmoothExponent = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Sp",1))
        settingsSmoothPoints = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Sn",3))
        settingsSmoothRadius = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Sr",numpy.inf))
        settingsSmoothFactor = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Sf",0.5))

        # Interpolate data on cell centers
        # TODO [Low]: Do not use subprocess here, use python directly
        exe = "foamUserToolInterpolate"
        arg = verbStr(verb) + " " \
            + "-Pp " + str(settingsInterpolateExponent) + " " \
            + "-Pn " + str(settingsInterpolatePoints) + " " \
            + "-Pr " + str(settingsInterpolateRadius) + " " \
            + settingsSmoothStr + " " \
            + "-Sp " + str(settingsSmoothExponent) + " " \
            + "-Sn " + str(settingsSmoothPoints) + " " \
            + "-Sr " + str(settingsSmoothRadius) + " " \
            + "-Sf " + str(settingsSmoothFactor) + " " \
            + "-i " + sourceComments + " " + tnodesComments + " " \
            + "-csd " + str(sourceColumns) + " " \
            + "-cso " + str(sourceOffset) + " " \
            + "-cno " + str(tnodesOffset) + " " \
            + sourcePath + " " + tnodesPath + " " + targetPath
        log = foamTimeLogPath + "/" + "foamUserToolInterpolate"
        subProc.run(exe, arg, log, verb)

        ## Loop over all data patches
        for p, patch in enumerate(patchList):

            # Read field information
            patchBaseType = patchTypeList[p][1]
            patchPrimitiveType, patchValued = foamClassInfo.patchBaseTypeToPrimitiveValued(patchBaseType)

            # Interpolate only data-valued patches
            if patchValued:

                # Print info
                if (verb > 0): print("interpolate data for field " + field + " on face centers of patch " + patch)

                # Read interpolation source information
                source = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("source",sub="intdata"))
                sourceOffset = foamDict.readInt(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("sourceOffset",sub="intdata"))
                sourceComments = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("sourceComments",sub="intdata"))
                sourcePath = foamTimeIntPath + "/" + source + ".dat"

                # Check interpolation source data for patch
                if not os.path.exists(sourcePath):
                   raise IOError("interpolation source data file \"" + sourcePath + "\" for patch " + patch + " of field " + field + " is missing")

                # Read interpolation target node information
                tnodes = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("tnodes",sub="intdata"))
                tnodesOffset = foamDict.readInt(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("tnodesOffset",sub="intdata"))
                tnodesComments = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("tnodesComments",sub="intdata"))
                tnodesPath = foamTimeRawPath + "/" + tnodes + "." + patch + ".dat"

                # Check interpolation target node data for patch
                if not os.path.exists(tnodesPath):
                   raise IOError("interpolation target node data file \"" + tnodesPath + "\" for patch " + patch + " of field " + field + " is missing")

                # Assemble target information
                target = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("target",sub="intdata"))
                targetPath = foamTimeRawPath + "/" + target + "." + patch + ".dat"

                # Parse interpolation data settings sub-dictionary for internal field
                intDataSettingsDict = foamDict.readSubDict(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("settings",sub="intdata"))

                # Read settings
                settingsInterpolateExponent = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Pp",3))
                settingsInterpolatePoints = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Pn",8))
                settingsInterpolateRadius = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Pr",numpy.inf))
                settingsSmooth = foamDict.readBool(foamDict.tryOrDefault(intDataSettingsDict, "S",False))
                settingsSmoothStr = ""
                if settingsSmooth: settingsSmoothStr = "-S"
                settingsSmoothExponent = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Sp",1))
                settingsSmoothPoints = foamDict.readInt(foamDict.tryOrDefault(intDataSettingsDict, "Sn",3))
                settingsSmoothRadius = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Sr",numpy.inf))
                settingsSmoothFactor = foamDict.readFloat(foamDict.tryOrDefault(intDataSettingsDict, "Sf",0.5))

                # Interpolate data on cell centers
                # TODO [Low]: Do not use subprocess here, use python directly
                exe = "foamUserToolInterpolate"
                arg = verbStr(verb) + " " \
                    + "-Pp " + str(settingsInterpolateExponent) + " " \
                    + "-Pn " + str(settingsInterpolatePoints) + " " \
                    + "-Pr " + str(settingsInterpolateRadius) + " " \
                    + settingsSmoothStr + " " \
                    + "-Sp " + str(settingsSmoothExponent) + " " \
                    + "-Sn " + str(settingsSmoothPoints) + " " \
                    + "-Sr " + str(settingsSmoothRadius) + " " \
                    + "-Sf " + str(settingsSmoothFactor) + " " \
                    + "-i " + sourceComments + " " + tnodesComments + " " \
                    + "-csd " + str(sourceColumns) + " " \
                    + "-cso " + str(sourceOffset) + " " \
                    + "-cno " + str(tnodesOffset) + " " \
                    + sourcePath + " " + tnodesPath + " " + targetPath
                log = foamTimeLogPath + "/" + "foamUserToolInterpolate"
                subProc.run(exe, arg, log, verb)


    # --- convert raw data to foam object format ---------------------------------- #

    # Create from raw if necessary
    if methodCreateFromRaw:

        # Print info
        if (verb > 0): print("create foam object from raw data of field " + field)

        # Assemble field info tuple
        fInfo = (field, fieldClass, fieldDimension)

        # Parse raw data dictionary for internal field
        rawDataDict = foamDict.readSubDict(foamDict.tryMandatory(fieldDictList[f], "rawdata"))

        # Read raw information and assemble raw info tuple
        raw = foamDict.readStr(foamDict.tryOrDefault(rawDataDict, "data",field))
        rawOffset = foamDict.readInt(foamDict.tryMandatory(rawDataDict, "offset"))
        rawComments = foamDict.readStr(foamDict.tryMandatory(rawDataDict, "comments"))
        rawPath = foamTimeRawPath + "/" + raw + "." + "internal" + ".dat"

        # Check raw data for internalField
        if not os.path.exists(rawPath):
            raise IOError("raw data file \"" + rawPath + "\" for field " + field + " is missing")

        # Assemble information for field conversion
        rInfo = (rawPath, rawOffset, rawComments)

        # Init patch info list
        pInfoList = []

        ## Loop over all data patches
        for p, patch in enumerate(patchList):

            # Read field information
            patchBaseType = patchTypeList[p][1]
            patchPrimitiveType, patchValued = foamClassInfo.patchBaseTypeToPrimitiveValued(patchBaseType)

            # Read raw information and assemble raw info tuple
            raw = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("data",sub="rawdata"))
            rawOffset = foamDict.readInt(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("offset",sub="rawdata"))
            rawComments = foamDict.readStr(foamFieldDictValIfPatchKeyNotPresentFromFieldUpdateDict("comments",sub="rawdata"))
            rawPath = foamTimeRawPath + "/" + raw + "." + patch + ".dat"

            # Check raw data for patch
            if patchValued and not os.path.exists(rawPath):
                raise IOError("raw data file \"" + rawPath + "\" for patch " + patch + " of field " + field + " is missing")

            # Assemble information for field conversion
            pData = (rawPath, rawOffset, rawComments)
            pInfo = (patch, patchBaseType, pData)
            pInfoList.append(pInfo)

        # Delete existing foam field files in time dir
        if (verb > 0): print("delete old files for field " + field)
        for path in [foamTimePath + "/" + field, foamTimePath + "/" + field + ".gz"]:
            if not os.path.isdir(path):
                if os.path.exists(path): os.remove(path)

        # Write raw data to foam field object
        fieldWriteFromRaw(time, fInfo, rInfo, pInfoList)



# --------------------------------------------------------------------------- #
# --- End of script --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Remove crap file
os.path.exists("PlyParser_FoamFileParser_parsetab.py") and os.remove("PlyParser_FoamFileParser_parsetab.py")

# Exit
sys.exit(0)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
