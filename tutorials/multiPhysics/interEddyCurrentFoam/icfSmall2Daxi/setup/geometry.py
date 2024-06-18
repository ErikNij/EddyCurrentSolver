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

v["all"]    = [0, 5, 35, 30]
v["fluid"]  = [12, 13, 25, 24]
v["vessel"] = [6, 8, 26, 25, 13, 12]
v["coila"]  = [15, 16, 22, 21]
v["free"]   = [0, 5, 35, 30, 24, 26, 8, 6]

s = dict()

for k in v.keys():

    s[k] = makeSketch(d, "inner", orient="xz")

    sketchPolyLine(s[k], v[k], blockMeshDict.v)

cs = par.coil_scale/par.geo_scale

cv    = dict()

cv[0] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                    -par.coil_bundle["z"]/2.0])
cv[1] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                    -par.coil_bundle["z"]/2.0])
cv[2] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                     par.coil_bundle["z"]/2.0])
cv[3] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                     par.coil_bundle["z"]/2.0])

s["coil"] = makeSketch(d, "coil", orient="xz",
                       base=(par.coils_origin[0], par.coils_origin[1],
                             par.coils_origin[2]))

sketchPolyLine(s["coil"], cv.keys(), cv)

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Bodies ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

bo = dict()

for k in s.keys():

    bo[k] = makeDoubleRevolveBody(k, s[k], angle=par.mesh_angle)

bo["conductor"] = makeFuseBody("conductor", [bo["fluid"],  bo["vessel"]])

bo["coils"] = makeOrthoArrayBody("coils", bo["coil"],
                                 (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

bo2D = dict()

for k in s.keys():

    bo2D[k] = makeDoubleExtrudeBody(k + "_2D", s[k], par.mesh_thickness)

bo2D["conductor"] = makeFuseBody("conductor_2D", [bo2D["fluid"], bo2D["vessel"]])

bo2D["coils"] = makeOrthoArrayBody("coils_2D", bo2D["coil"],
                                   (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

bo3D = dict()

for k in s.keys():

    bo3D[k] = makeRevolveBody(k + "_3D", s[k])

bo3D["conductor"] = makeFuseBody("conductor_3D", [bo3D["fluid"], bo3D["vessel"]])

bo3D["coils"] = makeOrthoArrayBody("coils_3D", bo3D["coil"],
                                   (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

shd = dict()

shd["front"] = (bo["all"], [4])
shd["back"] = (bo["all"], [7])
shd["infinity"] = (bo["all"], [1, 2, 3, 5, 6, 8])
shd["atmosphere"] = (bo["fluid"], [6, 8])
shd["sideWall"] = (bo["fluid"], [2, 5])
shd["bottomWall"] = (bo["fluid"], [1, 3])

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

    exportObj = bodyObj + bo2D.values() + bo3D.values() + shellObj

    #Export body and shell meshes
    exportMeshes(exportObj, par.dir_triSurface, __head__, scale=par.geo_scale)

    # Save document
    d.saveAs(par.dir_triSurface + "/" + __head__ + ".fcstd")

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
