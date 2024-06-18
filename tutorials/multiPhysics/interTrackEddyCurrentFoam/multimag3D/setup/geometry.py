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

s = dict()

s["inner"] = makeSketch(d, "inner", orient="xy", base=(0.0, 0.0, par.geo_z1))

sketchCircle(s["inner"], par.geo_r1)

s["outer"] = makeSketch(d, "outer", orient="xy", base=(0.0, 0.0, par.geo_z1))

sketchCircle(s["outer"], par.geo_r2)

cs = par.coil_scale/par.geo_scale

if par.coil_setup == "RMF":

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

    s["coil"] = makeSketch(d, "coil", orient="yz",
                           base=(par.coils_origin[0] + par.coils_step,
                                 par.coils_origin[1], par.coils_origin[2]))

    sketchPolyLine(s["coil"] , cvi.keys(), cvi,
                   fillet=(par.coil_path["r"] - par.coil_bundle["r"]/2.0))

    sketchPolyLine(s["coil"] , cvo.keys(), cvo,
                   fillet=(par.coil_path["r"] + par.coil_bundle["r"]/2.0))

elif par.coil_setup == "TMF":

    s["coil"] = makeSketch(d, "coil", orient="xy",
                           base=(par.coils_origin[0], par.coils_origin[1],
                                 par.coils_origin[2]))

    sketchCircle(s["coil"], cs*(par.coil_path["r"] - par.coil_bundle['r']/2.0))
    sketchCircle(s["coil"], cs*(par.coil_path["r"] + par.coil_bundle['r']/2.0))

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Bodies ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

bo = dict()

bo["fluid"] = makeExtrudeBody("fluid", s["inner"], par.geo_z3-par.geo_z1)
bo["above"] = makeExtrudeBody("above", s["outer"], par.geo_z4-par.geo_z1)
bo["below"] = makeExtrudeBody("below", s["outer"], par.geo_z0-par.geo_z1)

bo["buffer"] = makeCutBody("buffer", bo["above"], bo["fluid"])

bo["all"] = makeFuseBody("all", [bo["above"], bo["below"]])

bo["conductor"] = makeFuseBody("conductor", [bo["fluid"]])

if par.coil_setup == "RMF":

    bo["coil"] = makeDoubleExtrudeBody("coil", s["coil"],
                                       cs*par.coil_bundle["z"])

    bo["coils"] = makePolarArrayBody("coils", bo["coil"], par.coils_n)

elif par.coil_setup == "TMF":

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

shd["infinity"] = (bo["all"], [1, 2, 3, 4])
shd["fixedMesh"] = (bo["buffer"], [1, 2, 3])
shd["sideWall"] = (bo["fluid"], [1])
shd["bottomWall"] = (bo["fluid"], [2])
shd["trackedSurface"] = (bo["fluid"], [3])

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

    exportObj = bodyObj + shellObj

    #Export body abnd shell meshes
    exportMeshes(exportObj, par.dir_triSurface, __head__, scale=par.geo_scale)

    # Save document
    d.saveAs(par.dir_triSurface + "/" + __head__ + ".fcstd")

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
