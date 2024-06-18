#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# June 2016
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
geo_r0 =   1.0
geo_r1 =   2.0
geo_r2 =  10.0

a = 2.0**(-0.5)
b = 3.0**( 0.5)/3.0
f = 1.2

geo_a0 = a*geo_r0
geo_a1 = a*geo_r1
geo_a2 = a*geo_r2

geo_b0 = b*geo_r0
geo_b1 = b*geo_r1
geo_b2 = b*geo_r2

n_scale = 2.0

n0 = int(m.ceil(n_scale*(geo_r0)*10))
n1 = int(m.ceil(n_scale*(geo_r1-geo_r0)*10))
n2 = int(m.ceil(n_scale*((geo_r2-geo_r1)+(geo_r1-geo_r0))/6.0*10))

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set(  0, [        0,        0,        0])
d.vertices.set(  1, [   geo_r0,        0,        0])
d.vertices.set(  2, [ f*geo_a0, f*geo_a0,        0])
d.vertices.set(  3, [        0,   geo_r0,        0])
d.vertices.set(  4, [        0,        0,   geo_r0])
d.vertices.set(  5, [ f*geo_a0,        0, f*geo_a0])
d.vertices.set(  6, [ f*geo_b0, f*geo_b0, f*geo_b0])
d.vertices.set(  7, [        0, f*geo_a0, f*geo_a0])

d.vertices.set(  8, [   geo_r1,        0,        0])
d.vertices.set(  9, [   geo_a1,   geo_a1,        0])
d.vertices.set( 10, [        0,   geo_r1,        0])
d.vertices.set( 11, [   geo_a1,        0,   geo_a1])
d.vertices.set( 12, [   geo_b1,   geo_b1,   geo_b1])
d.vertices.set( 13, [        0,   geo_a1,   geo_a1])
d.vertices.set( 14, [        0,        0,   geo_r1])

d.vertices.set( 15, [   geo_r2,        0,        0])
d.vertices.set( 16, [   geo_a2,   geo_a2,        0])
d.vertices.set( 17, [        0,   geo_r2,        0])
d.vertices.set( 18, [   geo_a2,        0,   geo_a2])
d.vertices.set( 19, [   geo_b2,   geo_b2,   geo_b2])
d.vertices.set( 20, [        0,   geo_a2,   geo_a2])
d.vertices.set( 21, [        0,        0,   geo_r2])

# Blocks

d.blocks.set( 0, [  0,  1,  2,  3,  4,  5,  6,  7], zone="conductor")

d.blocks.set( 1, [  1,  8,  9,  2,  5, 11, 12,  6], zone="conductor")
d.blocks.set( 2, [  3,  2,  9, 10,  7,  6, 12, 13], zone="conductor")
d.blocks.set( 3, [  4,  5,  6,  7, 14, 11, 12, 13], zone="conductor")

d.blocks.set( 4, [  8, 15, 16,  9, 11, 18, 19, 12], zone="background")
d.blocks.set( 5, [ 10,  9, 16, 17, 13, 12, 19, 20], zone="background")
d.blocks.set( 6, [ 14, 11, 12, 13, 21, 18, 19, 20], zone="background")

# Distributions

d.blocks.distribution.set(0, "x", n0)
d.blocks.distribution.set(0, "y", n0)
d.blocks.distribution.set(0, "z", n0)

d.blocks.distribution.set(1, "x", n1)

d.blocks.distribution.set(4, "x", n2)

# Gradings

d.blocks.grading.set( 4, [8.0, 1.0, 1.0])

# Boundary faces

d.boundaryFaces.set( 0, "mirror_x", 0, "x-")
d.boundaryFaces.set( 1, "mirror_x", 2, "x-")
d.boundaryFaces.set( 2, "mirror_x", 3, "x-")
d.boundaryFaces.set( 3, "mirror_x", 5, "x-")
d.boundaryFaces.set( 4, "mirror_y", 6, "x-")

d.boundaryFaces.set( 5, "mirror_y", 0, "y-")
d.boundaryFaces.set( 6, "mirror_y", 1, "y-")
d.boundaryFaces.set( 7, "mirror_y", 3, "y-")
d.boundaryFaces.set( 8, "mirror_y", 4, "y-")
d.boundaryFaces.set( 9, "mirror_y", 6, "y-")

d.boundaryFaces.set(10, "mirror_z", 0, "z-")
d.boundaryFaces.set(11, "mirror_z", 1, "z-")
d.boundaryFaces.set(12, "mirror_z", 2, "z-")
d.boundaryFaces.set(13, "mirror_z", 4, "z-")
d.boundaryFaces.set(14, "mirror_z", 5, "z-")

d.boundaryFaces.set(15, "infinity", 4, "x+")
d.boundaryFaces.set(16, "infinity", 5, "y+")
d.boundaryFaces.set(17, "infinity", 6, "z+")

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

    d.arc(  8,  9, 10)
    d.arc( 10,  9,  8)

    d.arc(  8, 11, 14)
    d.arc( 14, 11,  8)

    d.arc( 10, 13, 14)
    d.arc( 14, 13, 10)

    d.arc(  9, 12, 14)
    d.arc( 11, 12, 10)
    d.arc( 13, 12,  8)

    d.arc( 15, 16, 17)
    d.arc( 17, 16, 15)

    d.arc( 15, 18, 21)
    d.arc( 21, 18, 15)

    d.arc( 17, 20, 21)
    d.arc( 21, 20, 17)

    d.arc( 16, 19, 21)
    d.arc( 18, 19, 17)
    d.arc( 20, 19, 15)

if d.subDict("boundary"):

    if d.boundarySubDict("mirror_x", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_y", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_z", "patch"):

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
