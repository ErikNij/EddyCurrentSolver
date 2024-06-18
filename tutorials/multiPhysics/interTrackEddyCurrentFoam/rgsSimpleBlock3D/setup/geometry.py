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

from foamTools.freecad import (makeSketch, sketchCircle, sketchPolyLine,
                               makeFuseBody, makeCutBody,
                               makeMirrorBody, makeExtrudeBody, makeRevolveBody,
                               makeDoubleExtrudeBody, makeDoubleRevolveBody,
                               makeOrthoArrayBody, makePolarArrayBody,
                               makeFaceShell, exportMeshes)

import FreeCAD

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par
import blockMeshDict

FreeCAD.newDocument(__head__)
FreeCAD.setActiveDocument(__head__)
d = FreeCAD.activeDocument()

# --------------------------------------------------------------------------- #
# --- Sketches -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

v = dict()

v["all"]         = [30, 34, 49, 45]
v["above"]       = [0, 4, 49, 45]
v["below"]       = [30, 34, 4, 0]
v["conductor0"]  = [1, 2, 17, 16]
v["conductor1"]  = [2, 3, 18, 17]

t = dict()

t["all"]         = par.mesh_X["inf"]
t["above"]       = par.mesh_X["inf"]
t["below"]       = par.mesh_X["inf"]
t["conductor0"]  = par.geo_X["C"]
t["conductor1"]  = par.geo_X["C"]

s = dict()

for k in v.keys():

    s[k] = makeSketch(d, "inner", orient="yz")

    vrt = dict(blockMeshDict.v)

    # Transform points into xy
    for l in vrt: vrt[l] = [vrt[l][1], vrt[l][2], vrt[l][0]]

    sketchPolyLine(s[k], v[k], vrt)

#cs = par.coil_scale/par.geo_scale

#cv    = dict()

#cv[0] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                    #-par.coil_bundle["z"]/2.0])
#cv[1] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                    #-par.coil_bundle["z"]/2.0])
#cv[2] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                     #par.coil_bundle["z"]/2.0])
#cv[3] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                     #par.coil_bundle["z"]/2.0])

#s["coil"] = makeSketch(d, "coil", orient="xz",
                       #base=(par.coils_origin[0], par.coils_origin[1],
                             #par.coils_origin[2]))

#sketchPolyLine(s["coil"], cv.keys(), cv)

cs = par.coil_scale/par.geo_scale

cvi    = dict()
cvo    = dict()

cvi[0] = cs*np.array([-par.coil_path["x"] + par.coil_bundle["r"]/2.0,
                      -par.coil_path["y"] + par.coil_bundle["r"]/2.0])
cvi[1] = cs*np.array([ par.coil_path["x"] - par.coil_bundle["r"]/2.0,
                      -par.coil_path["y"] + par.coil_bundle["r"]/2.0])
cvi[2] = cs*np.array([ par.coil_path["x"] - par.coil_bundle["r"]/2.0,
                       par.coil_path["y"] - par.coil_bundle["r"]/2.0])
cvi[3] = cs*np.array([-par.coil_path["x"] + par.coil_bundle["r"]/2.0,
                       par.coil_path["y"] - par.coil_bundle["r"]/2.0])

cvo[0] = cs*np.array([-par.coil_path["x"] - par.coil_bundle["r"]/2.0,
                      -par.coil_path["y"] - par.coil_bundle["r"]/2.0])
cvo[1] = cs*np.array([ par.coil_path["x"] + par.coil_bundle["r"]/2.0,
                      -par.coil_path["y"] - par.coil_bundle["r"]/2.0])
cvo[2] = cs*np.array([ par.coil_path["x"] + par.coil_bundle["r"]/2.0,
                       par.coil_path["y"] + par.coil_bundle["r"]/2.0])
cvo[3] = cs*np.array([-par.coil_path["x"] - par.coil_bundle["r"]/2.0,
                       par.coil_path["y"] + par.coil_bundle["r"]/2.0])

s["coil"] = makeSketch(d, "coil", orient="xy",
                       base=(par.coils_origin[0], par.coils_origin[1],
                             par.coils_origin[2]))

sketchPolyLine(s["coil"] , cvi.keys(), cvi,
               fillet=(par.coil_path["r"] - par.coil_bundle["r"]/2.0))

sketchPolyLine(s["coil"] , cvo.keys(), cvo,
               fillet=(par.coil_path["r"] + par.coil_bundle["r"]/2.0))

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Bodies ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

bo = dict()

for k in s.keys():

    if k is not "coil": bo[k] = makeDoubleExtrudeBody(k, s[k], 2.0 * t[k])

bo["conductor"] = makeFuseBody("conductor", [bo["conductor0"],
                                             bo["conductor1"]])

bo["fluid"] = makeFuseBody("fluid", [bo["conductor"]])

bo["buffer"] = makeCutBody("buffer", bo["above"], bo["fluid"])

bo["coil"] = makeDoubleExtrudeBody("coil", s["coil"],
                                   cs*par.coil_bundle["z"])

bo["coils"] = makeOrthoArrayBody("coils", bo["coil"],
                                 (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

shd = dict()

shd["infinity"] = (bo["all"], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
shd["fixedMesh"] = (bo["buffer"], [1, 2, 3, 4, 5, 6, 11, 12, 13, 18])
shd["sideWall"] = (bo["fluid"], [1, 3, 4, 8, 10, 12, 14, 16])
shd["bottomWall"] = (bo["fluid"], [5, 9, 11, 15])
shd["trackedSurface"] = (bo["fluid"], [2, 6, 7, 13])

sh = dict()

for k in shd.keys():

    sh[k] = makeFaceShell(d, k, shd[k])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

bodyObj = bo.values()
shellObj = sh.values()

def main():

    exportObj = bodyObj + bo.values() + shellObj

    #Export body and shell meshes
    exportMeshes(exportObj, par.dir_triSurface, __head__, scale=par.geo_scale)

    # Save document
    d.saveAs(par.dir_triSurface + "/" + __head__ + ".fcstd")

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
