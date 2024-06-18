#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.blockMeshDict import *
import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale = 1e-3

# Background mesh
geo_x1 = -150.0
geo_x2 =  150.0

geo_dy =    1.25 # 5 mm / (2*2)
geo_y1 = geo_dy/2.0
geo_y2 = geo_y1 + geo_dy

geo_z1 =  -60.0
geo_z2 =  140.0

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

d.blocks.set( 0, [0, 1, 5, 4, 2, 3, 7, 6], [60, 1, 40], zone="background")

d.boundaryFaces.set(1, "front", 0, "y-")
d.boundaryFaces.set(2, "back", 0, "y+")
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

    if d.boundarySubDict("front", "empty"):

        d.boundaryFaces.write()

    if d.boundarySubDict("back", "empty"):

        d.boundaryFaces.write()

    if d.boundarySubDict("infinity", "patch"):

        d.boundaryFaces.write()

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
