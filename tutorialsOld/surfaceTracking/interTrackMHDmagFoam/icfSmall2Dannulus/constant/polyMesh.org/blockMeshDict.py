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

geo_xSlip = 1.0
geo_x0 =  85.0
geo_x1 = geo_x0 +  1.0
geo_x2 = geo_x1 + 14.0
geo_x3 = geo_x2 + 15.0
geo_x4 = 300.0

geo_z0 = geo_x1 - geo_x0
geo_z1 = 200.0
geo_z2 = 235.0
geo_z3 = 383.0
geo_z4 = geo_z1 + geo_z3

axi_phi        = 5.0
axi_phiRad     = m.pi * axi_phi/180.0
axi_phiHalf    = axi_phiRad/2.0
axi_phiHalfSin = m.sin(axi_phiHalf)
axi_phiHalfCos = m.cos(axi_phiHalf)

axi_xSlip = geo_xSlip * axi_phiHalfCos
axi_x0 = geo_x0 * axi_phiHalfCos
axi_x1 = geo_x1 * axi_phiHalfCos
axi_x2 = geo_x2 * axi_phiHalfCos
axi_x3 = geo_x3 * axi_phiHalfCos
axi_x4 = geo_x4 * axi_phiHalfCos

axi_ySlip = geo_xSlip * axi_phiHalfSin
axi_y0 = geo_x0 * axi_phiHalfSin
axi_y1 = geo_x1 * axi_phiHalfSin
axi_y2 = geo_x2 * axi_phiHalfSin
axi_y3 = geo_x3 * axi_phiHalfSin
axi_y4 = geo_x4 * axi_phiHalfSin

n_scale = 0.25

n_x0 = int(m.ceil(n_scale*geo_x0))
n_x1 = int(m.ceil(n_scale*(geo_x1-geo_x0)))
n_x2 = int(m.ceil(n_scale*(geo_x2-geo_x1)))
n_x3 = int(m.ceil(n_scale*(geo_x3-geo_x2)))
n_x4 = int(m.ceil(n_scale*(geo_x4-geo_x3)))

n_z0 = int(m.ceil(n_scale*geo_z0))
n_z1 = int(m.ceil(n_scale*(geo_z1-geo_z0)))
n_z2 = int(m.ceil(n_scale*geo_z2))
n_z3 = int(m.ceil(n_scale*(geo_z3-geo_z2)))
n_z4 = int(m.ceil(n_scale*(geo_z4-geo_z3)))

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set( 0, [    0.0,     0.0, -geo_z1])
d.vertices.set( 1, [    0.0,     0.0, -geo_z0])
d.vertices.set( 2, [    0.0,     0.0,     0.0])
d.vertices.set( 3, [    0.0,     0.0,  geo_z2])
d.vertices.set( 4, [    0.0,     0.0,  geo_z3])
d.vertices.set( 5, [    0.0,     0.0,  geo_z4])

d.vertices.set( 6, [ axi_x0, -axi_y0, -geo_z1])
d.vertices.set( 7, [ axi_x0, -axi_y0, -geo_z0])
d.vertices.set( 8, [ axi_x0, -axi_y0,     0.0])
d.vertices.set( 9, [ axi_x0, -axi_y0,  geo_z2])
d.vertices.set(10, [ axi_x0, -axi_y0,  geo_z3])
d.vertices.set(11, [ axi_x0, -axi_y0,  geo_z4])

d.vertices.set(12, [ axi_x1, -axi_y1, -geo_z1])
d.vertices.set(13, [ axi_x1, -axi_y1, -geo_z0])
d.vertices.set(14, [ axi_x1, -axi_y1,     0.0])
d.vertices.set(15, [ axi_x1, -axi_y1,  geo_z2])
d.vertices.set(16, [ axi_x1, -axi_y1,  geo_z3])
d.vertices.set(17, [ axi_x1, -axi_y1,  geo_z4])

d.vertices.set(18, [ axi_x2, -axi_y2, -geo_z1])
d.vertices.set(19, [ axi_x2, -axi_y2, -geo_z0])
d.vertices.set(20, [ axi_x2, -axi_y2,     0.0])
d.vertices.set(21, [ axi_x2, -axi_y2,  geo_z2])
d.vertices.set(22, [ axi_x2, -axi_y2,  geo_z3])
d.vertices.set(23, [ axi_x2, -axi_y2,  geo_z4])

d.vertices.set(24, [ axi_x3, -axi_y3, -geo_z1])
d.vertices.set(25, [ axi_x3, -axi_y3, -geo_z0])
d.vertices.set(26, [ axi_x3, -axi_y3,     0.0])
d.vertices.set(27, [ axi_x3, -axi_y3,  geo_z2])
d.vertices.set(28, [ axi_x3, -axi_y3,  geo_z3])
d.vertices.set(29, [ axi_x3, -axi_y3,  geo_z4])

d.vertices.set(30, [ axi_x4, -axi_y4, -geo_z1])
d.vertices.set(31, [ axi_x4, -axi_y4, -geo_z0])
d.vertices.set(32, [ axi_x4, -axi_y4,     0.0])
d.vertices.set(33, [ axi_x4, -axi_y4,  geo_z2])
d.vertices.set(34, [ axi_x4, -axi_y4,  geo_z3])
d.vertices.set(35, [ axi_x4, -axi_y4,  geo_z4])

d.vertices.set(36, [ axi_x0,  axi_y0, -geo_z1])
d.vertices.set(37, [ axi_x0,  axi_y0, -geo_z0])
d.vertices.set(38, [ axi_x0,  axi_y0,     0.0])
d.vertices.set(39, [ axi_x0,  axi_y0,  geo_z2])
d.vertices.set(40, [ axi_x0,  axi_y0,  geo_z3])
d.vertices.set(41, [ axi_x0,  axi_y0,  geo_z4])

d.vertices.set(42, [ axi_x1,  axi_y1, -geo_z1])
d.vertices.set(43, [ axi_x1,  axi_y1, -geo_z0])
d.vertices.set(44, [ axi_x1,  axi_y1,     0.0])
d.vertices.set(45, [ axi_x1,  axi_y1,  geo_z2])
d.vertices.set(46, [ axi_x1,  axi_y1,  geo_z3])
d.vertices.set(47, [ axi_x1,  axi_y1,  geo_z4])

d.vertices.set(48, [ axi_x2,  axi_y2, -geo_z1])
d.vertices.set(49, [ axi_x2,  axi_y2, -geo_z0])
d.vertices.set(50, [ axi_x2,  axi_y2,     0.0])
d.vertices.set(51, [ axi_x2,  axi_y2,  geo_z2])
d.vertices.set(52, [ axi_x2,  axi_y2,  geo_z3])
d.vertices.set(53, [ axi_x2,  axi_y2,  geo_z4])

d.vertices.set(54, [ axi_x3,  axi_y3, -geo_z1])
d.vertices.set(55, [ axi_x3,  axi_y3, -geo_z0])
d.vertices.set(56, [ axi_x3,  axi_y3,     0.0])
d.vertices.set(57, [ axi_x3,  axi_y3,  geo_z2])
d.vertices.set(58, [ axi_x3,  axi_y3,  geo_z3])
d.vertices.set(59, [ axi_x3,  axi_y3,  geo_z4])

d.vertices.set(60, [ axi_x4,  axi_y4, -geo_z1])
d.vertices.set(61, [ axi_x4,  axi_y4, -geo_z0])
d.vertices.set(62, [ axi_x4,  axi_y4,     0.0])
d.vertices.set(63, [ axi_x4,  axi_y4,  geo_z2])
d.vertices.set(64, [ axi_x4,  axi_y4,  geo_z3])
d.vertices.set(65, [ axi_x4,  axi_y4,  geo_z4])

d.vertices.set( 0, [ axi_xSlip, -axi_ySlip, -geo_z1])
d.vertices.set( 1, [ axi_xSlip, -axi_ySlip, -geo_z0])
d.vertices.set( 2, [ axi_xSlip, -axi_ySlip,     0.0])
d.vertices.set( 3, [ axi_xSlip, -axi_ySlip,  geo_z2])
d.vertices.set( 4, [ axi_xSlip, -axi_ySlip,  geo_z3])
d.vertices.set( 5, [ axi_xSlip, -axi_ySlip,  geo_z4])

d.vertices.set(66, [ axi_xSlip,  axi_ySlip, -geo_z1])
d.vertices.set(67, [ axi_xSlip,  axi_ySlip, -geo_z0])
d.vertices.set(68, [ axi_xSlip,  axi_ySlip,     0.0])
d.vertices.set(69, [ axi_xSlip,  axi_ySlip,  geo_z2])
d.vertices.set(70, [ axi_xSlip,  axi_ySlip,  geo_z3])
d.vertices.set(71, [ axi_xSlip,  axi_ySlip,  geo_z4])

#d.blocks.set( 0, [  0,  6, 36,  0,  1,  7, 37,  1], zone="region_static")
#d.blocks.set( 1, [  1,  7, 37,  1,  2,  8, 38,  2], zone="region_static")
#d.blocks.set( 2, [  2,  8, 38,  2,  3,  9, 39,  3], zone="region_fluid")
#d.blocks.set( 3, [  3,  9, 39,  3,  4, 10, 40,  4], zone="region_dynamic")
#d.blocks.set( 4, [  4, 10, 40,  4,  5, 11, 41,  5], zone="region_static")

d.blocks.set( 5, [  6, 12, 42, 36,  7, 13, 43, 37], zone="region_static")
d.blocks.set( 6, [  7, 13, 43, 37,  8, 14, 44, 38], zone="region_static")
d.blocks.set( 7, [  8, 14, 44, 38,  9, 15, 45, 39], zone="region_dynamic")
d.blocks.set( 8, [  9, 15, 45, 39, 10, 16, 46, 40], zone="region_dynamic")
d.blocks.set( 9, [ 10, 16, 46, 40, 11, 17, 47, 41], zone="region_static")

d.blocks.set(10, [ 12, 18, 48, 42, 13, 19, 49, 43], zone="region_static")
d.blocks.set(11, [ 13, 19, 49, 43, 14, 20, 50, 44], zone="region_static")
d.blocks.set(12, [ 14, 20, 50, 44, 15, 21, 51, 45], zone="region_static")
d.blocks.set(13, [ 15, 21, 51, 45, 16, 22, 52, 46], zone="region_static")
d.blocks.set(14, [ 16, 22, 52, 46, 17, 23, 53, 47], zone="region_static")

d.blocks.set(15, [ 18, 24, 54, 48, 19, 25, 55, 49], zone="region_static")
d.blocks.set(16, [ 19, 25, 55, 49, 20, 26, 56, 50], zone="region_static")
d.blocks.set(17, [ 20, 26, 56, 50, 21, 27, 57, 51], zone="region_static")
d.blocks.set(18, [ 21, 27, 57, 51, 22, 28, 58, 52], zone="region_static")
d.blocks.set(19, [ 22, 28, 58, 52, 23, 29, 59, 53], zone="region_static")

d.blocks.set(20, [ 24, 30, 60, 54, 25, 31, 61, 55], zone="region_static")
d.blocks.set(21, [ 25, 31, 61, 55, 26, 32, 62, 56], zone="region_static")
d.blocks.set(22, [ 26, 32, 62, 56, 27, 33, 63, 57], zone="region_static")
d.blocks.set(23, [ 27, 33, 63, 57, 28, 34, 64, 58], zone="region_static")
d.blocks.set(24, [ 28, 34, 64, 58, 29, 35, 65, 59], zone="region_static")

d.blocks.set( 0, [  0,  6, 36, 66,  1,  7, 37, 67], zone="region_static")
d.blocks.set( 1, [  1,  7, 37, 67,  2,  8, 38, 68], zone="region_static")
d.blocks.set( 2, [  2,  8, 38, 68,  3,  9, 39, 69], zone="region_fluid")
d.blocks.set( 3, [  3,  9, 39, 69,  4, 10, 40, 70], zone="region_dynamic")
d.blocks.set( 4, [  4, 10, 40, 70,  5, 11, 41, 71], zone="region_static")

d.blocks.distribution.set( 0, "x", n_x0)
d.blocks.distribution.set( 5, "x", n_x1)
d.blocks.distribution.set(10, "x", n_x2)
d.blocks.distribution.set(15, "x" ,n_x3)
d.blocks.distribution.set(20, "x", n_x4)

d.blocks.distribution.set( 0, "z", n_z1)
d.blocks.distribution.set( 1, "z", n_z0)
d.blocks.distribution.set( 2, "z", n_z2)
d.blocks.distribution.set( 3, "z", n_z3)
d.blocks.distribution.set( 4, "z", n_z4)

d.boundaryFaces.set( 0, "axis", [  0,  1,  2,  3,  4], "x-")
d.boundaryFaces.set( 1, "front", [  0,  1,  2,  3,  4, \
                                    5,  6,  7,  8,  9, \
                                   10, 11, 12, 13, 14, \
                                   15, 16, 17, 18, 19, \
                                   20, 21, 22, 23, 24], "y-")
d.boundaryFaces.set( 2, "back",  [  0,  1,  2,  3,  4, \
                                    5,  6,  7,  8,  9, \
                                   10, 11, 12, 13, 14, \
                                   15, 16, 17, 18, 19, \
                                   20, 21, 22, 23, 24], "y+")
d.boundaryFaces.set( 3, "infinity", [  0,  5, 10, 15, 20], "z-")
d.boundaryFaces.set( 4, "infinity", [  4,  9, 14, 19, 24], "z+")
d.boundaryFaces.set( 5, "infinity", [ 20, 21, 22, 23, 24], "x+")

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

    if d.boundarySubDict("axis", "symmetryPlane"):

        d.boundaryFaces.write()

    if d.boundarySubDict("front", "wedge"):

        d.boundaryFaces.write()

    if d.boundarySubDict("back", "wedge"):

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
