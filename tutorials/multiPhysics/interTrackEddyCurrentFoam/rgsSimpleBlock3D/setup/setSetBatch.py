#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
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

sys.path.append(os.environ["FOAM_USER_TOOLS"] + "/" + "python")
sys.path.append("/usr/lib/freecad/lib")

import math as m
import numpy as np

from foamTools.setSetBatch import setSetBatch

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par
import geometry

maxLength = 60
tol = 2e-4

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    ssb = setSetBatch(maxLength=maxLength)

    ssb.rename(par.dir_system + "/" + "setSetBatch")

    if ssb.group("All"):

        ssb.cellSet("all", "all")

    if ssb.group("Bodies"):

        for e in geometry.bodyObj:

            name = str(e.Label)

            relDirPath = os.path.relpath(par.dir_triSurface)

            relFilePath = relDirPath + "/" + "geometry_" + name + ".stl"

            ssb.cellSet(name, "surfaceToTopo", relFilePath, tol=tol, inside=True)

    if ssb.group("Shells"):

        for e in geometry.shellObj:

            name = str(e.Label)

            relDirPath = os.path.relpath(par.dir_triSurface)

            relFilePath = relDirPath + "/" + "geometry_" + name + ".stl"

            ssb.faceSet(name, "surfaceToTopo", relFilePath, tol=tol)

    if ssb.group("Regions"):

        ssb.cellSet("region_fluid", "setToTopo",
                    add=["body_fluid"])

        #ssb.cellSet("region_buffer", "setToTopo",
                    #add=["body_buffer"])
        ssb.cellSet("region_buffer", "setToTopo",
                    add=["body_above"],
                    delete=["body_fluid"])

        ssb.cellSet("region_conductor", "setToTopo",
                    add=["body_fluid"])

        ssb.faceSet("regions", "cellBoundaryToTopo",
                    add=["region_fluid", "region_buffer", "region_conductor"],
                    delete=["all"])

    if ssb.group("Patches"):

        ssb.faceSet("patch_infinity", "setToTopo", add=["shell_infinity"])
        ssb.faceSet("patch_fixedMesh", "setToTopo", add=["shell_fixedMesh"])
        ssb.faceSet("patch_bottomWall", "setToTopo", add=["shell_bottomWall"])
        ssb.faceSet("patch_sideWall", "setToTopo", add=["shell_sideWall"])
        ssb.faceSet("patch_trackedSurface", "setToTopo", add=["shell_trackedSurface"])

    if ssb.group("Materials"):

        ssb.cellSet("material_conductor0", "setToTopo", add=["body_conductor0"])
        ssb.cellSet("material_conductor1", "setToTopo", add=["body_conductor1"])

    ssb.quit()

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
