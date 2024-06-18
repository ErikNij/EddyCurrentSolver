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

from foamTools.blockMeshDict import expansion_de_e, expansion_n_ds, blockMeshDict

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

a = 2.0**(-0.5)
f = par.mesh_f

d = blockMeshDict(mesh=par.mesh)

# --------------------------------------------------------------------------- #
# --- Vertices -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.vertices.set(  0, [            0.0,              0,     par.geo_z0])
d.vertices.set(  1, [     par.geo_r0,              0,     par.geo_z0])
d.vertices.set(  2, [ f*a*par.geo_r0, f*a*par.geo_r0,     par.geo_z0])
d.vertices.set(  3, [            0.0,     par.geo_r0,     par.geo_z0])
d.vertices.set(  4, [     par.geo_r1,              0,     par.geo_z0])
d.vertices.set(  5, [   a*par.geo_r1,   a*par.geo_r1,     par.geo_z0])
d.vertices.set(  6, [            0.0,     par.geo_r1,     par.geo_z0])
d.vertices.set(  7, [     par.geo_r2,              0,     par.geo_z0])
d.vertices.set(  8, [   a*par.geo_r2,   a*par.geo_r2,     par.geo_z0])
d.vertices.set(  9, [            0.0,     par.geo_r2,     par.geo_z0])

# --------------------------------------------------------------------------- #

baseVertices = [ i for i in range(len(d.vertices.labels)) ]

d.vertices.copyTranslate( 10, baseVertices, [0.0, 0.0, par.geo_z1-par.geo_z0])
d.vertices.copyTranslate( 20, baseVertices, [0.0, 0.0, par.geo_z2-par.geo_z0])
d.vertices.copyTranslate( 30, baseVertices, [0.0, 0.0, par.geo_z3-par.geo_z0])
d.vertices.copyTranslate( 40, baseVertices, [0.0, 0.0, par.geo_z4-par.geo_z0])

# --------------------------------------------------------------------------- #
# --- Blocks ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.set(  0, [  0,  1,  2,  3, 10, 11, 12, 13], zone="core")
d.blocks.set(  1, [  1,  4,  5,  2, 11, 14, 15, 12], zone="core")
d.blocks.set(  2, [  3,  2,  5,  6, 13, 12, 15, 16], zone="core")

d.blocks.set(  3, [  4,  7,  8,  5, 14, 17, 18, 15], zone="annulus")
d.blocks.set(  4, [  6,  5,  8,  9, 16, 15, 18, 19], zone="annulus")

# --------------------------------------------------------------------------- #

baseBlocks = [ i for i in range(len(d.blocks.labels)) ]

d.blocks.copyShiftVerticeLabels(  10, baseBlocks,  10)
d.blocks.copyShiftVerticeLabels(  20, baseBlocks,  20)
d.blocks.copyShiftVerticeLabels(  30, baseBlocks,  30)

# --------------------------------------------------------------------------- #
# --- Expansion ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

lr0 = par.geo_r0
lr1 = par.geo_r1 - par.geo_r0
lr2 = par.geo_r2 - par.geo_r1

lz0 = par.geo_z1-par.geo_z0
lz1 = par.geo_z2-par.geo_z1
lz2 = par.geo_z3-par.geo_z2
lz3 = par.geo_z4-par.geo_z3

nr0 = int(m.ceil(par.mesh_scale*lr0*4))
nr1 = int(m.ceil(par.mesh_scale*lr1*5))
#nr2 = int(m.ceil(par.mesh_scale*lr2*2))

#nz0 = int(m.ceil(par.mesh_scale*lz0*2))
nz1 = int(m.ceil(par.mesh_scale*lz1*5))
nz2 = int(m.ceil(par.mesh_scale*lz2*5))
#nz3 = int(m.ceil(par.mesh_scale*lz3*2))

e1 = 0.5
e2ds = expansion_de_e(nr1, e1, lr1)
e2de = 9.0*e2ds*lr2/90.0
e2 = e2de/e2ds
nr2 = expansion_n_ds(e2, e2ds, lr2)

nz0 = nr2
nz3 = nr2

# --------------------------------------------------------------------------- #
# --- Distribution ---------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.distribution.set(  0, "x", nr0)
d.blocks.distribution.set(  0, "y", nr0)
d.blocks.distribution.set(  1, "x", nr1)
d.blocks.distribution.set(  3, "x", nr2)

d.blocks.distribution.set(  0, "z", nz0)
d.blocks.distribution.set( 10, "z", nz1)
d.blocks.distribution.set( 20, "z", nz2)
d.blocks.distribution.set( 30, "z", nz3)

# --------------------------------------------------------------------------- #
# --- Grading --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.blocks.grading.set( 1, [ e1, 1.0,    1.0])
d.blocks.grading.set( 3, [ e2, 1.0,    1.0])

d.blocks.grading.set( 0, [1.0, 1.0, 1.0/e2])
d.blocks.grading.set(10, [1.0, 1.0, 1.0/e1])
d.blocks.grading.set(20, [1.0, 1.0,     e1])
d.blocks.grading.set(30, [1.0, 1.0,     e2])

# --------------------------------------------------------------------------- #
# --- Boundary -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.boundaryFaces.set(  0, "mirror_x",  0, "x-")
d.boundaryFaces.set(  1, "mirror_x",  2, "x-")
d.boundaryFaces.set(  2, "mirror_x",  4, "x-")
d.boundaryFaces.set( 10, "mirror_x", 10, "x-")
d.boundaryFaces.set( 11, "mirror_x", 12, "x-")
d.boundaryFaces.set( 12, "mirror_x", 14, "x-")
d.boundaryFaces.set( 20, "mirror_x", 20, "x-")
d.boundaryFaces.set( 21, "mirror_x", 22, "x-")
d.boundaryFaces.set( 22, "mirror_x", 24, "x-")
d.boundaryFaces.set( 30, "mirror_x", 30, "x-")
d.boundaryFaces.set( 31, "mirror_x", 32, "x-")
d.boundaryFaces.set( 32, "mirror_x", 34, "x-")

d.boundaryFaces.set(  3, "mirror_y",  0, "y-")
d.boundaryFaces.set(  4, "mirror_y",  1, "y-")
d.boundaryFaces.set(  6, "mirror_y",  3, "y-")
d.boundaryFaces.set( 13, "mirror_y", 10, "y-")
d.boundaryFaces.set( 14, "mirror_y", 11, "y-")
d.boundaryFaces.set( 16, "mirror_y", 13, "y-")
d.boundaryFaces.set( 23, "mirror_y", 20, "y-")
d.boundaryFaces.set( 24, "mirror_y", 21, "y-")
d.boundaryFaces.set( 26, "mirror_y", 23, "y-")
d.boundaryFaces.set( 33, "mirror_y", 30, "y-")
d.boundaryFaces.set( 34, "mirror_y", 31, "y-")
d.boundaryFaces.set( 36, "mirror_y", 33, "y-")

d.boundaryFaces.set( 100, "infinity",  0, "z-")
d.boundaryFaces.set( 101, "infinity",  1, "z-")
d.boundaryFaces.set( 102, "infinity",  2, "z-")
d.boundaryFaces.set( 103, "infinity",  3, "z-")
d.boundaryFaces.set( 104, "infinity",  4, "z-")

d.boundaryFaces.set( 130, "infinity", 30, "z+")
d.boundaryFaces.set( 131, "infinity", 31, "z+")
d.boundaryFaces.set( 132, "infinity", 32, "z+")
d.boundaryFaces.set( 133, "infinity", 33, "z+")
d.boundaryFaces.set( 134, "infinity", 34, "z+")

d.boundaryFaces.set( 200, "infinity",  3, "x+")
d.boundaryFaces.set( 201, "infinity",  4, "y+")
d.boundaryFaces.set( 210, "infinity", 13, "x+")
d.boundaryFaces.set( 211, "infinity", 14, "y+")
d.boundaryFaces.set( 220, "infinity", 23, "x+")
d.boundaryFaces.set( 221, "infinity", 24, "y+")
d.boundaryFaces.set( 230, "infinity", 33, "x+")
d.boundaryFaces.set( 231, "infinity", 34, "y+")

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

        d.arc(  4,  5,  6)
        d.arc(  6,  5,  4)

        d.arc( 14, 15, 16)
        d.arc( 16, 15, 14)

        d.arc( 24, 25, 26)
        d.arc( 26, 25, 24)

        d.arc( 34, 35, 36)
        d.arc( 36, 35, 34)

        d.arc( 44, 45, 46)
        d.arc( 46, 45, 44)

        d.arc(  7,  8,  9)
        d.arc(  9,  8,  7)

        d.arc( 17, 18, 19)
        d.arc( 19, 18, 17)

        d.arc( 27, 28, 29)
        d.arc( 29, 28, 27)

        d.arc( 37, 38, 39)
        d.arc( 39, 38, 37)

        d.arc( 47, 48, 49)
        d.arc( 49, 48, 47)

        pass

    if d.subDict("boundary"):

        if d.boundarySubDict("mirror_x", "patch"):

            d.boundaryFaces.write()

        if d.boundarySubDict("mirror_y", "patch"):

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
