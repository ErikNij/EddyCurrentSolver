#!/usr/bin/python
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

geo_scale = 1

# Background mesh
geo_x1 =  0.0
geo_x2 =  1.0

geo_y1 =  0.0
geo_y2 =  1.0

geo_z1 =  0.0
geo_z2 =  1.0

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

d.blocks.set( 0, [0, 1, 5, 4, 2, 3, 7, 6], [20, 20, 20])

d.boundaryFaces.set(1, "bottomWall", 0, "y-")
d.boundaryFaces.set(2, "trackedSurface", 0, "y+")
d.boundaryFaces.set(3, "sideWalls", 0, "x-")
d.boundaryFaces.set(4, "sideWalls", 0, "x+")
d.boundaryFaces.set(5, "sideWalls", 0, "z-")
d.boundaryFaces.set(6, "sideWalls", 0, "z+")

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

    if d.boundarySubDict("bottomWall", "wall"):

        d.boundaryFaces.write()

    if d.boundarySubDict("trackedSurface", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("sideWalls", "wall"):

        d.boundaryFaces.write()

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

