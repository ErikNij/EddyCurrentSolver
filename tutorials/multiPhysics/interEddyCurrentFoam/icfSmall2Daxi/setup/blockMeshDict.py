#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (par.beckstein@hzdr.de)

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

from foamTools.blockMeshDict import blockMeshDict

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

d = blockMeshDict(mesh=par.mesh)

# --------------------------------------------------------------------------- #
# --- Vertices -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def vmg(v): return np.linalg.norm(v)
def vr0(r): return np.array([r, 0.0])
def v0z(z): return np.array([0.0, z])
def vrz(r, z): return np.array([r, z])
def vsa(s, a): return s * np.array([m.cos(a), m.sin(a)])
def vwg(al, wl=None):
    if wl is None: wl = np.ones(len(al))
    a = np.zeros(2)
    for i, ai in enumerate(al): a = np.add(a, wl[i]*ai)
    return a/np.sum(wl)

v      = dict()

v[  0] = v0z(par.mesh_Z["C"] - par.mesh_Z["inf"])
v[  1] = vrz(par.geo_R[0], v[0][1])
v[  2] = vrz(par.geo_R[1], v[0][1])
v[  3] = vrz(par.geo_R[2], v[0][1])
v[  4] = vrz(par.geo_R[3], v[0][1])
v[  5] = vrz(par.mesh_R["inf"], v[0][1])

v[  6] = v0z(par.geo_Z[0])
v[  7] = vrz(par.geo_R[0], v[6][1])
v[  8] = vrz(par.geo_R[1], v[6][1])
v[  9] = vrz(par.geo_R[2], v[6][1])
v[ 10] = vrz(par.geo_R[3], v[6][1])
v[ 11] = vrz(par.mesh_R["inf"], v[6][1])

v[ 12] = v0z(par.geo_Z[1])
v[ 13] = vrz(par.geo_R[0], v[12][1])
v[ 14] = vrz(par.geo_R[1], v[12][1])
v[ 15] = vrz(par.geo_R[2], v[12][1])
v[ 16] = vrz(par.geo_R[3], v[12][1])
v[ 17] = vrz(par.mesh_R["inf"], v[12][1])

v[ 18] = v0z(par.geo_Z[2])
v[ 19] = vrz(par.geo_R[0], v[18][1])
v[ 20] = vrz(par.geo_R[1], v[18][1])
v[ 21] = vrz(par.geo_R[2], v[18][1])
v[ 22] = vrz(par.geo_R[3], v[18][1])
v[ 23] = vrz(par.mesh_R["inf"], v[18][1])

v[ 24] = v0z(par.geo_Z[3])
v[ 25] = vrz(par.geo_R[0], v[24][1])
v[ 26] = vrz(par.geo_R[1], v[24][1])
v[ 27] = vrz(par.geo_R[2], v[24][1])
v[ 28] = vrz(par.geo_R[3], v[24][1])
v[ 29] = vrz(par.mesh_R["inf"], v[24][1])

v[ 30] = v0z(par.mesh_Z["C"] + par.mesh_Z["inf"])
v[ 31] = vrz(par.geo_R[0], v[30][1])
v[ 32] = vrz(par.geo_R[1], v[30][1])
v[ 33] = vrz(par.geo_R[2], v[30][1])
v[ 34] = vrz(par.geo_R[3], v[30][1])
v[ 35] = vrz(par.mesh_R["inf"], v[30][1])

# --------------------------------------------------------------------------- #

# Copy vertice dict and its containing vertices
V = v.copy()
for i in V: V[i] = v[i].copy()

# Axis point shift based on inner point p
def dv0z(v, p): r = par.mesh_R["axis"]; return np.array([r, r*(p[1]-v[1])/p[0]])

# Shift axis points
V[ 0] += dv0z(V[ 0], V[ 1])
V[ 6] += dv0z(V[ 6], V[ 7])
V[12] += dv0z(V[12], V[13])
V[18] += dv0z(V[18], V[19])
V[24] += dv0z(V[24], V[25])
V[30] += dv0z(V[30], V[31])

# --------------------------------------------------------------------------- #

for i in V: d.vertices.set(i, V[i])

# --------------------------------------------------------------------------- #
# --- Blocks ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.set(  0, [  0,  1,  7,  6], zone="free_external")
d.blocks.set(  1, [  1,  2,  8,  7], zone="free_external")
d.blocks.set(  2, [  2,  3,  9,  8], zone="free_external")
d.blocks.set(  3, [  3,  4, 10,  9], zone="free_external")
d.blocks.set(  4, [  4,  5, 11, 10], zone="free_external")

d.blocks.set(  5, [  6,  7, 13, 12], zone="vessel")
d.blocks.set(  6, [  7,  8, 14, 13], zone="vessel")
d.blocks.set(  7, [  8,  9, 15, 14], zone="free_external")
d.blocks.set(  8, [  9, 10, 16, 15], zone="free_external")
d.blocks.set(  9, [ 10, 11, 17, 16], zone="free_external")

d.blocks.set( 10, [ 12, 13, 19, 18], zone="fluid")
d.blocks.set( 11, [ 13, 14, 20, 19], zone="vessel")
d.blocks.set( 12, [ 14, 15, 21, 20], zone="free_external")
d.blocks.set( 13, [ 15, 16, 22, 21], zone="coila")
d.blocks.set( 14, [ 16, 17, 23, 22], zone="free_external")

d.blocks.set( 15, [ 18, 19, 25, 24], zone="fluid")
d.blocks.set( 16, [ 19, 20, 26, 25], zone="vessel")
d.blocks.set( 17, [ 20, 21, 27, 26], zone="free_external")
d.blocks.set( 18, [ 21, 22, 28, 27], zone="free_external")
d.blocks.set( 19, [ 22, 23, 29, 28], zone="free_external")

d.blocks.set( 20, [ 24, 25, 31, 30], zone="free_external")
d.blocks.set( 21, [ 25, 26, 32, 31], zone="free_external")
d.blocks.set( 22, [ 26, 27, 33, 32], zone="free_external")
d.blocks.set( 23, [ 27, 28, 34, 33], zone="free_external")
d.blocks.set( 24, [ 28, 29, 35, 34], zone="free_external")

# --------------------------------------------------------------------------- #
# --- Expansion ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def n(l): s = 0.1; return int(m.ceil(abs(s*par.mesh_scale*l)))

nr                  = dict()
nr["free_side"]     = n(v[35][0] - v[34][0])
nr["coila"]         = n(v[34][0] - v[33][0])
nr["free_gap"]      = n(v[33][0] - v[32][0])
nr["vessel"]        = n(v[32][0] - v[31][0])
nr["fluid"]         = n(v[31][0] - v[30][0])

nz                  = dict()
nz["free_top"]      = n(v[30][1] - v[24][1])
nz["fluid_air"]     = n(v[24][1] - v[18][1])
nz["fluid_melt"]    = n(v[18][1] - v[12][1])
nz["free_bottom"]   = n(v[0][1] - v[6][1])

# --------------------------------------------------------------------------- #
# --- Distribution ---------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.distribution.set( 24, "x", nr["free_side"])
d.blocks.distribution.set( 23, "x", nr["coila"])
d.blocks.distribution.set( 22, "x", nr["free_gap"])
d.blocks.distribution.set( 21, "x", nr["vessel"])
d.blocks.distribution.set( 20, "x", nr["fluid"])

d.blocks.distribution.set( 20, "y", nz["free_top"])
d.blocks.distribution.set( 15, "y", nz["fluid_air"])
d.blocks.distribution.set( 10, "y", nz["fluid_melt"])
d.blocks.distribution.set(  0, "y", nz["free_bottom"])

# --------------------------------------------------------------------------- #
# --- Grading --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#d.blocks.grading.set( 0, [1/3.0, 1.0, 1.0])

# --------------------------------------------------------------------------- #
# --- Boundary -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.boundaryFaces.set(  0, "axis", [  0,  5, 10, 15, 20], "x-")
d.boundaryFaces.set(  1, "front", d.blocks.labels, "z+")
d.boundaryFaces.set(  2, "back", d.blocks.labels, "z-")
d.boundaryFaces.set(  3, "infinity", [  4,  9, 14, 19, 24], "x+")
d.boundaryFaces.set(  4, "infinity", [ 20, 21, 22, 23, 24], "y+")
d.boundaryFaces.set(  5, "infinity", [  0,  1,  2,  3,  4], "y-")

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

        if d.boundarySubDict("axis", "patch"):

            d.boundaryFaces.write()

        if d.boundarySubDict("front", "empty"):

            d.boundaryFaces.write()

        if d.boundarySubDict("back", "empty"):

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
