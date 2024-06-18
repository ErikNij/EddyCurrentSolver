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

geo_scale = 1e-2

# Background mesh
geo_rS = 0.01
geo_r0 = 0.625
geo_r1 = 1.25
geo_r2 = 5.0

a = 2.0**(-0.5)
f = 1.2

axi_phi        = 5.0
axi_phiRad     = m.pi * axi_phi/180.0
axi_phiHalf    = axi_phiRad/2.0
axi_phiHalfSin = m.sin(axi_phiHalf)
axi_phiHalfCos = m.cos(axi_phiHalf)

axi_xS = geo_rS * axi_phiHalfCos
axi_x0 = geo_r0 * axi_phiHalfCos
axi_x1 = geo_r1 * axi_phiHalfCos
axi_x2 = geo_r2 * axi_phiHalfCos

axi_yS = geo_rS * axi_phiHalfSin
axi_y0 = geo_r0 * axi_phiHalfSin
axi_y1 = geo_r1 * axi_phiHalfSin
axi_y2 = geo_r2 * axi_phiHalfSin

n_scale = 4.0

n0 = int(m.ceil(n_scale*(geo_r0)*5))
n1 = int(m.ceil(n_scale*(geo_r1-geo_r0)*10))
n2 = int(m.ceil(n_scale*((geo_r2-geo_r1)+(geo_r1-geo_r0))/6.0*10))

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set(  0, [     axi_xS,     axi_yS,          0])
d.vertices.set(  1, [     axi_x0,     axi_y0,          0])
d.vertices.set(  2, [ f*a*axi_x0, f*a*axi_y0, f*a*geo_r0])
d.vertices.set(  3, [     axi_xS,     axi_yS,     geo_r0])
d.vertices.set(  4, [     axi_xS,    -axi_yS,          0])
d.vertices.set(  5, [     axi_x0,    -axi_y0,          0])
d.vertices.set(  6, [ f*a*axi_x0,-f*a*axi_y0, f*a*geo_r0])
d.vertices.set(  7, [     axi_xS,    -axi_yS,     geo_r0])

d.vertices.set(  8, [     axi_x1,     axi_y1,          0])
d.vertices.set(  9, [   a*axi_x1,   a*axi_y1,   a*geo_r1])
d.vertices.set( 10, [     axi_xS,     axi_yS,     geo_r1])
d.vertices.set( 11, [     axi_x1,    -axi_y1,          0])
d.vertices.set( 12, [   a*axi_x1,  -a*axi_y1,   a*geo_r1])
d.vertices.set( 13, [     axi_xS,    -axi_yS,     geo_r1])

d.vertices.set( 14, [     axi_x2,     axi_y2,          0])
d.vertices.set( 15, [   a*axi_x2,   a*axi_y2,   a*geo_r2])
d.vertices.set( 16, [     axi_xS,     axi_yS,     geo_r2])
d.vertices.set( 17, [     axi_x2,    -axi_y2,          0])
d.vertices.set( 18, [   a*axi_x2,  -a*axi_y2,   a*geo_r2])
d.vertices.set( 19, [     axi_xS,    -axi_yS,     geo_r2])

# Blocks

d.blocks.set( 0, [  0,  1,  2,  3,  4,  5,  6,  7], zone="conductor")

d.blocks.set( 1, [  1,  8,  9,  2,  5, 11, 12,  6], zone="conductor")
d.blocks.set( 2, [  3,  2,  9, 10,  7,  6, 12, 13], zone="conductor")

d.blocks.set( 3, [  8, 14, 15,  9, 11, 17, 18, 12], zone="background")
d.blocks.set( 4, [ 10,  9, 15, 16, 13, 12, 18, 19], zone="background")

# Distributions

d.blocks.distribution.set(0, "x", n0)
d.blocks.distribution.set(0, "y", n0)
d.blocks.distribution.set(0, "z", 1)

d.blocks.distribution.set(1, "x", n1)

d.blocks.distribution.set(3, "x", n2)

# Gradings

d.blocks.grading.set( 1, [0.25, 1.0, 1.0])
d.blocks.grading.set( 3, [20.0, 1.0, 1.0])

# Boundary faces

d.boundaryFaces.set(  0, "axis", 0, "x-")
d.boundaryFaces.set(  1, "axis", 2, "x-")
d.boundaryFaces.set(  2, "axis", 4, "x-")

d.boundaryFaces.set(  3, "front", 0, "z-")
d.boundaryFaces.set(  4, "front", 1, "z-")
d.boundaryFaces.set(  5, "front", 2, "z-")
d.boundaryFaces.set(  6, "front", 3, "z-")
d.boundaryFaces.set(  7, "front", 4, "z-")

d.boundaryFaces.set(  8, "back", 0, "z+")
d.boundaryFaces.set(  9, "back", 1, "z+")
d.boundaryFaces.set( 10, "back", 2, "z+")
d.boundaryFaces.set( 11, "back", 3, "z+")
d.boundaryFaces.set( 12, "back", 4, "z+")

d.boundaryFaces.set( 13, "infinity", 3, "x+")
d.boundaryFaces.set( 14, "infinity", 4, "y+")

d.boundaryFaces.set( 15, "mirror_z", 0, "y-")
d.boundaryFaces.set( 16, "mirror_z", 1, "y-")
d.boundaryFaces.set( 17, "mirror_z", 3, "y-")

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

    d.arc( 11, 12, 13)
    d.arc( 13, 12, 11)

    d.arc( 14, 15, 16)
    d.arc( 16, 15, 14)

    d.arc( 17, 18, 19)
    d.arc( 19, 18, 17)

if d.subDict("boundary"):

    if d.boundarySubDict("axis", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("front", "wedge"):

        d.boundaryFaces.write()

    if d.boundarySubDict("back", "wedge"):

        d.boundaryFaces.write()

    if d.boundarySubDict("infinity", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_z", "patch"):

        d.boundaryFaces.write()
    pass

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
