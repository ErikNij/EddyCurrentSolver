#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.blockMeshDict import blockMeshDict
import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale = 1e-3

# Background mesh
geo_x1 = -180.0
geo_x2 =  180.0

geo_y1 = -220.0
geo_y2 =  220.0

geo_z1 = -115.0 # d1 =  35.0 # d2 =  35.0
geo_z2 =  160.0 # d1 = 115.0 # d1 =  70.0

n_scale = 0.1

n_x = int(m.ceil(n_scale*(geo_x2-geo_x1)))
n_y = int(m.ceil(n_scale*(geo_y2-geo_y1)))
n_z = int(m.ceil(n_scale*(geo_z2-geo_z1)))

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set(  0, [geo_x1, geo_y1, geo_z1])
d.vertices.set(  1, [geo_x2, geo_y1, geo_z1])
d.vertices.set(  2, [geo_x1, geo_y1, geo_z2])
d.vertices.set(  3, [geo_x2, geo_y1, geo_z2])
d.vertices.set(  4, [geo_x1, geo_y2, geo_z1])
d.vertices.set(  5, [geo_x2, geo_y2, geo_z1])
d.vertices.set(  6, [geo_x1, geo_y2, geo_z2])
d.vertices.set(  7, [geo_x2, geo_y2, geo_z2])

d.blocks.set( 0, [0, 1, 5, 4, 2, 3, 7, 6], [n_x, n_y, n_z], zone="background")

d.boundaryFaces.set(1, "infinity", 0, "y-")
d.boundaryFaces.set(2, "infinity", 0, "y+")
d.boundaryFaces.set(3, "infinity", 0, "x-")
d.boundaryFaces.set(4, "infinity", 0, "x+")
d.boundaryFaces.set(5, "infinity", 0, "z-")
d.boundaryFaces.set(6, "infinity", 0, "z+")

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.header(geo_scale)

# --------------------------------------------------------------------------- #

if d.subDict("vertices"):

    d.vertices.write()

if d.subDict("blocks"):

    d.blocks.write()

if d.subDict("edges"):

    pass

if d.subDict("boundary"):

    if d.boundarySubDict("infinity", "patch"):

        d.boundaryFaces.write()

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
