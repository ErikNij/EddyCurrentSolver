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

import math as m
import numpy as np

from foamTools.blockMeshDict import expansion_n_ds, blockMeshDict

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

d = blockMeshDict(mesh=par.mesh)

# --------------------------------------------------------------------------- #
# --- Vertices -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def vxy0(x, y): return np.array([x, y, 0.0])

v     = dict()

v[ 0] = vxy0(0.0, par.mesh_Y["inf"])
v[ 1] = vxy0(0.0, par.geo_Y["C"])
v[ 2] = vxy0(0.0, par.geo_Y["shift"])
v[ 3] = vxy0(0.0, -par.geo_Y["C"])
v[ 4] = vxy0(0.0, -par.mesh_Y["inf"])

v[ 5] = vxy0(par.geo_X["C"], par.mesh_Y["inf"])
v[ 6] = vxy0(par.geo_X["C"], par.geo_Y["C"])
v[ 7] = vxy0(par.geo_X["C"], par.geo_Y["shift"])
v[ 8] = vxy0(par.geo_X["C"], -par.geo_Y["C"])
v[ 9] = vxy0(par.geo_X["C"], -par.mesh_Y["inf"])

v[10] = vxy0(par.mesh_X["inf"], par.mesh_Y["inf"])
v[11] = vxy0(par.mesh_X["inf"], par.geo_Y["C"])
v[12] = vxy0(par.mesh_X["inf"], 0.0)
v[13] = vxy0(par.mesh_X["inf"], -par.geo_Y["C"])
v[14] = vxy0(par.mesh_X["inf"], -par.mesh_Y["inf"])

# --------------------------------------------------------------------------- #

shift = len(v)

for l in range(0,15):

    v[shift+l] = np.array(v[l])
    v[shift+l][2] = par.geo_Z["C"]

for l in [2, 7, 12]: v[shift+l][1] = -par.geo_Y["shift"]
for l in [12]: v[shift+l][1] = 0.0

shift = len(v)

for l in range(0,15):

    v[shift+l] = np.array(v[l])
    v[shift+l][2] = -par.mesh_Z["inf-"]

for l in [2, 7, 12]: v[shift+l][1] = 0.0

shift = len(v)

for l in range(0,15):

    v[shift+l] = np.array(v[l])
    v[shift+l][2] = par.mesh_Z["inf+"]

for l in [2, 7, 12]: v[shift+l][1] = 0.0

# --------------------------------------------------------------------------- #

for i in v: d.vertices.set(i, v[i])

# --------------------------------------------------------------------------- #
# --- Blocks ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.set(  0, [  0,  1,  6,  5, 15, 16, 21 , 20], zone="free")
d.blocks.set(  1, [  1,  2,  7,  6, 16, 17, 22 , 21], zone="conductor")
d.blocks.set(  2, [  2,  3,  8,  7, 17, 18, 23 , 22], zone="conductor")
d.blocks.set(  3, [  3,  4,  9,  8, 18, 19, 24 , 23], zone="free")

d.blocks.set(  4, [  5,  6, 11, 10, 20, 21, 26 , 25], zone="free")
d.blocks.set(  5, [  6,  7, 12, 11, 21, 22, 27 , 26], zone="free")
d.blocks.set(  6, [  7,  8, 13, 12, 22, 23, 28 , 27], zone="free")
d.blocks.set(  7, [  8,  9, 14, 13, 23, 24, 29 , 28], zone="free")

# --------------------------------------------------------------------------- #

baseBlocks = [ i for i in range(len(d.blocks.labels)) ]

d.blocks.copyShiftVerticeLabels( 8, baseBlocks, [ 30, 30, 30, 30,-15,-15,-15,-15])
d.blocks.copyShiftVerticeLabels(16, baseBlocks, [ 15, 15, 15, 15, 30, 30, 30, 30])

# --------------------------------------------------------------------------- #
# --- Expansion ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

ex                 = dict()
ex["free"]         = 8.0

ey                 = dict()
ey["free_y+"]      = ex["free"]/(v[10][0] - v[5][0]) * (v[0][1] - v[1][1])
ey["free_y-"]      = ex["free"]/(v[10][0] - v[5][0]) * (v[3][1] - v[4][1])

ez                 = dict()
ez["free_bottom"]  = ex["free"]/(v[10][0] - v[5][0]) * (v[0][2] - v[30][2])
ez["free_top"]     = ex["free"]/(v[10][0] - v[5][0]) * (v[45][2] - v[15][2])

# --------------------------------------------------------------------------- #

def n(l): s = 0.2; return int(m.ceil(abs(s*par.mesh_scale*l)))

nx                 = dict()
nx["conductor"]    = n(v[5][0] - v[0][0])
nx["free"]         = expansion_n_ds(ex["free"],
                                    (v[5][0] - v[0][0])/nx["conductor"],
                                    v[10][0] - v[5][0])

ny                 = dict()
ny["free_y+"]      = expansion_n_ds(ey["free_y+"],
                                    (v[5][0] - v[0][0])/nx["conductor"],
                                    v[0][1] - v[1][1])
ny["conductor_y+"] = n(v[1][1] - 0.0)
ny["conductor_y-"] = n(0.0 - v[3][1])
ny["free_y-"]      = expansion_n_ds(ey["free_y-"],
                                    (v[5][0] - v[0][0])/nx["conductor"],
                                    v[3][1] - v[4][1])

nz                 = dict()
nz["conductor"]    = n(v[15][2] - v[0][2])
nz["free_bottom"]  = expansion_n_ds(ez["free_bottom"],
                                    (v[5][0] - v[0][0])/nx["conductor"],
                                    v[0][2] - v[30][2])
nz["free_top"]     = expansion_n_ds(ez["free_top"],
                                    (v[5][0] - v[0][0])/nx["conductor"],
                                    v[45][2] - v[15][2])

# --------------------------------------------------------------------------- #
# --- Distribution ---------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.distribution.set( 0, "y", nx["conductor"])
d.blocks.distribution.set( 4, "y", nx["free"])

d.blocks.distribution.set( 0, "x", ny["free_y+"])
d.blocks.distribution.set( 1, "x", ny["conductor_y+"])
d.blocks.distribution.set( 2, "x", ny["conductor_y-"])
d.blocks.distribution.set( 3, "x", ny["free_y-"])

d.blocks.distribution.set( 0, "z", nz["conductor"])
d.blocks.distribution.set( 8, "z", nz["free_bottom"])
d.blocks.distribution.set(16, "z", nz["free_top"])

# --------------------------------------------------------------------------- #
# --- Grading --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.grading.set( 4, "y", ex["free"])

d.blocks.grading.set( 0, "x", 1.0/ey["free_y+"])
d.blocks.grading.set( 3, "x", ey["free_y+"])

d.blocks.grading.set( 8, "z", 1.0/ez["free_bottom"])
d.blocks.grading.set(16, "z", ez["free_bottom"])

# --------------------------------------------------------------------------- #
# --- Boundary -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.boundaryFaces.set( 0, "mirror_x", range( 0, 4), "y-")
d.boundaryFaces.set( 1, "mirror_x", range( 8,12), "y-")
d.boundaryFaces.set( 2, "mirror_x", range(16,20), "y-")

d.boundaryFaces.set( 3, "infinity", range( 8,16), "z-")
d.boundaryFaces.set( 4, "infinity", range(16,24), "z+")
d.boundaryFaces.set( 5, "infinity", [  0,  4,  8, 12, 16, 20], "x-")
d.boundaryFaces.set( 6, "infinity", [  3,  7, 11, 15, 19, 23], "x+")
d.boundaryFaces.set( 7, "infinity", range( 4, 8), "y+")
d.boundaryFaces.set( 8, "infinity", range(12,16), "y+")
d.boundaryFaces.set( 9, "infinity", range(20,24), "y+")

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    d.rename(par.dir_polyMesh + "/" + "blockMeshDict")

    d.header(par.geo_scale)

    if d.subDict("vertices"):

        d.vertices.write()

    if d.subDict("blocks"):

        d.blocks.write()

    if d.subDict("edges"):

        pass

    if d.subDict("boundary"):

        if d.boundarySubDict("mirror_x", "patch"):

            d.boundaryFaces.write()

        if d.boundarySubDict("infinity", "patch"):

            d.boundaryFaces.write()

        pass

    if d.subDict("mergePatchPairs"):

        pass

    d.footer()

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
