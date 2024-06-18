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
geo_x3 =  2.0
geo_x4 =  3.0

geo_y1 =  0.0
geo_y2 =  1.0

geo_dz =  0.0223607 * 2.0
geo_z1 = geo_dz/2.0
geo_z2 = geo_z1 + geo_dz

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set(  0, [geo_x1, geo_y1, geo_z1])
d.vertices.set(  1, [geo_x2, geo_y1, geo_z1])
d.vertices.set(  2, [geo_x1, geo_y1, geo_z2])
d.vertices.set(  3, [geo_x2, geo_y1, geo_z2])

d.vertices.set(  4, [geo_x1, geo_y2, geo_z1])
d.vertices.set(  5, [geo_x3, geo_y2, geo_z1])
d.vertices.set(  6, [geo_x1, geo_y2, geo_z2])
d.vertices.set(  7, [geo_x3, geo_y2, geo_z2])

d.vertices.set(  8, [geo_x4, geo_y1, geo_z1])
d.vertices.set(  9, [geo_x4, geo_y2, geo_z1])
d.vertices.set( 10, [geo_x4, geo_y1, geo_z2])
d.vertices.set( 11, [geo_x4, geo_y2, geo_z2])

d.blocks.set( 0, [0, 1, 5, 4, 2, 3, 7, 6], [20, 20, 1])
d.blocks.set( 1, [1, 8, 9, 5, 3, 10, 11, 7], [20, 20, 1])

d.boundaryFaces.set(1, "bottom", [0, 1], "y-")
d.boundaryFaces.set(2, "top", [0, 1], "y+")
d.boundaryFaces.set(3, "left", 0, "x-")
d.boundaryFaces.set(4, "right", 1, "x+")
d.boundaryFaces.set(5, "frontAndBack", [0, 1], "z-")
d.boundaryFaces.set(6, "frontAndBack", [0, 1], "z+")

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

    if d.boundarySubDict("bottom", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("top", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("left", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("right", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("frontAndBack", "empty"):

        d.boundaryFaces.write()

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

