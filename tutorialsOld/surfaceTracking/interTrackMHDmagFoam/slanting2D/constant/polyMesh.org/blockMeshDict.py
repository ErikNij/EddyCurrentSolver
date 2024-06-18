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

geo_lx0 =  35.0
geo_lx1 = 145.0
geo_lx2 =  45.0

#geo_ly0 =  10.0
geo_ly0 =   0.625

geo_lz0 =  60.0
geo_lz1 = 140.0
geo_lz2 =  50.0

geo_fz0 =  10.0
geo_fz1 =  30.0
geo_fz2 =  20.0

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")



d.vertices.set( 0, [-geo_lx1, -geo_ly0, -geo_lz0])
d.vertices.set( 1, [-geo_lx0, -geo_ly0, -geo_lz0])
d.vertices.set( 2, [ geo_lx0, -geo_ly0, -geo_lz0])
d.vertices.set( 3, [ geo_lx1, -geo_ly0, -geo_lz0])
d.vertices.set( 4, [-geo_lx1,  geo_ly0, -geo_lz0])
d.vertices.set( 5, [-geo_lx0,  geo_ly0, -geo_lz0])
d.vertices.set( 6, [ geo_lx0,  geo_ly0, -geo_lz0])
d.vertices.set( 7, [ geo_lx1,  geo_ly0, -geo_lz0])

d.vertices.set( 8, [-geo_lx1, -geo_ly0,      0.0])
d.vertices.set( 9, [-geo_lx0, -geo_ly0,      0.0])
d.vertices.set(10, [ geo_lx0, -geo_ly0,      0.0])
d.vertices.set(11, [ geo_lx1, -geo_ly0,      0.0])
d.vertices.set(12, [-geo_lx1,  geo_ly0,      0.0])
d.vertices.set(13, [-geo_lx0,  geo_ly0,      0.0])
d.vertices.set(14, [ geo_lx0,  geo_ly0,      0.0])
d.vertices.set(15, [ geo_lx1,  geo_ly0,      0.0])

d.vertices.set(16, [-geo_lx1, -geo_ly0,  geo_fz0])
d.vertices.set(17, [-geo_lx0, -geo_ly0,  geo_fz0])
d.vertices.set(18, [ geo_lx0, -geo_ly0,  geo_fz1])
d.vertices.set(19, [ geo_lx1, -geo_ly0,  geo_fz1])
d.vertices.set(20, [-geo_lx1,  geo_ly0,  geo_fz0])
d.vertices.set(21, [-geo_lx0,  geo_ly0,  geo_fz0])
d.vertices.set(22, [ geo_lx0,  geo_ly0,  geo_fz1])
d.vertices.set(23, [ geo_lx1,  geo_ly0,  geo_fz1])

d.vertices.set(24, [-geo_lx1, -geo_ly0,  geo_lz1])
d.vertices.set(25, [-geo_lx0, -geo_ly0,  geo_lz1])
d.vertices.set(26, [ geo_lx0, -geo_ly0,  geo_lz1])
d.vertices.set(27, [ geo_lx1, -geo_ly0,  geo_lz1])
d.vertices.set(28, [-geo_lx1,  geo_ly0,  geo_lz1])
d.vertices.set(29, [-geo_lx0,  geo_ly0,  geo_lz1])
d.vertices.set(30, [ geo_lx0,  geo_ly0,  geo_lz1])
d.vertices.set(31, [ geo_lx1,  geo_ly0,  geo_lz1])

d.blocks.set(1, [  0,  1,  5,  4,  8,  9, 13, 12], zone="region_dynamic")
d.blocks.set(2, [  1,  2,  6,  5,  9, 10, 14, 13], zone="region_dynamic")
d.blocks.set(3, [  2,  3,  7,  6, 10, 11, 15, 14], zone="region_dynamic")
d.blocks.set(4, [  8,  9, 13, 12, 16, 17, 21, 20], zone="region_dynamic")
d.blocks.set(5, [  9, 10, 14, 13, 17, 18, 22, 21], zone="region_fluid")
d.blocks.set(6, [ 10, 11, 15, 14, 18, 19, 23, 22], zone="region_dynamic")
d.blocks.set(7, [ 16, 17, 21, 20, 24, 25, 29, 28], zone="region_dynamic")
d.blocks.set(8, [ 17, 18, 22, 21, 25, 26, 30, 29], zone="region_dynamic")
d.blocks.set(9, [ 18, 19, 23, 22, 26, 27, 31, 30], zone="region_dynamic")

d.blocks.distribution.set(5, [35, 1, 10])

d.blocks.distribution.set(4, "x", 55)
d.blocks.distribution.set(6, "x", 55)

d.blocks.distribution.set(2, "z", 30)
d.blocks.distribution.set(8, "z", 60)

d.boundaryFaces.set(1, "front", [1, 2, 3, 4, 5, 6, 7, 8, 9], "y-")
d.boundaryFaces.set(2, "back", [1, 2, 3, 4, 5, 6, 7, 8, 9], "y+")
d.boundaryFaces.set(3, "infinity", [1, 2, 3], "z-")
d.boundaryFaces.set(4, "infinity", [7, 8, 9], "z+")
d.boundaryFaces.set(5, "infinity", [1, 4, 7], "x-")
d.boundaryFaces.set(6, "infinity", [3, 6, 9], "x+")



d.vertices.set(16, [-geo_lx1, -geo_ly0,  geo_fz2])
d.vertices.set(20, [-geo_lx1,  geo_ly0,  geo_fz2])

d.vertices.set(19, [ geo_lx1, -geo_ly0,  geo_fz2])
d.vertices.set(23, [ geo_lx1,  geo_ly0,  geo_fz2])

d.vertices.set(32, [-geo_lx2, -geo_ly0, -geo_lz0])
d.vertices.set(33, [-geo_lx2,  geo_ly0, -geo_lz0])
d.vertices.set(34, [-geo_lx2, -geo_ly0,      0.0])
d.vertices.set(35, [-geo_lx2,  geo_ly0,      0.0])
d.vertices.set(36, [-geo_lx2, -geo_ly0,  geo_fz2])
d.vertices.set(37, [-geo_lx2,  geo_ly0,  geo_fz2])
d.vertices.set(38, [-geo_lx2, -geo_ly0,  geo_lz2])
d.vertices.set(39, [-geo_lx2,  geo_ly0,  geo_lz2])
d.vertices.set(40, [-geo_lx2, -geo_ly0,  geo_lz1])
d.vertices.set(41, [-geo_lx2,  geo_ly0,  geo_lz1])

d.vertices.set(42, [ geo_lx2, -geo_ly0, -geo_lz0])
d.vertices.set(43, [ geo_lx2,  geo_ly0, -geo_lz0])
d.vertices.set(44, [ geo_lx2, -geo_ly0,      0.0])
d.vertices.set(45, [ geo_lx2,  geo_ly0,      0.0])
d.vertices.set(46, [ geo_lx2, -geo_ly0,  geo_fz2])
d.vertices.set(47, [ geo_lx2,  geo_ly0,  geo_fz2])
d.vertices.set(48, [ geo_lx2, -geo_ly0,  geo_lz2])
d.vertices.set(49, [ geo_lx2,  geo_ly0,  geo_lz2])
d.vertices.set(50, [ geo_lx2, -geo_ly0,  geo_lz1])
d.vertices.set(51, [ geo_lx2,  geo_ly0,  geo_lz1])

d.vertices.set(52, [-geo_lx1, -geo_ly0,  geo_lz2])
d.vertices.set(53, [-geo_lx1,  geo_ly0,  geo_lz2])
d.vertices.set(54, [-geo_lx0, -geo_ly0,  geo_lz2])
d.vertices.set(55, [-geo_lx0,  geo_ly0,  geo_lz2])
d.vertices.set(56, [ geo_lx0, -geo_ly0,  geo_lz2])
d.vertices.set(57, [ geo_lx0,  geo_ly0,  geo_lz2])
d.vertices.set(58, [ geo_lx1, -geo_ly0,  geo_lz2])
d.vertices.set(59, [ geo_lx1,  geo_ly0,  geo_lz2])

d.blocks.set( 1, [  0, 32, 33,  4,  8, 34, 35, 12], zone="region_static")
d.blocks.set( 2, [ 32,  1,  5, 33, 34,  9, 13, 35], zone="region_static")
d.blocks.set( 3, [  1,  2,  6,  5,  9, 10, 14, 13], zone="region_static")
d.blocks.set( 4, [  2, 42, 43,  6, 10, 44, 45, 14], zone="region_static")
d.blocks.set( 5, [ 42,  3,  7, 43, 44, 11, 15, 45], zone="region_static")
d.blocks.set( 6, [  8, 34, 35, 12, 16, 36, 37, 20], zone="region_static")
d.blocks.set( 7, [ 34,  9, 13, 35, 36, 17, 21, 37], zone="region_dynamic")
d.blocks.set( 8, [  9, 10, 14, 13, 17, 18, 22, 21], zone="region_fluid")
d.blocks.set( 9, [ 10, 44, 45, 14, 18, 46, 47, 22], zone="region_dynamic")
d.blocks.set(10, [ 44, 11, 15, 45, 46, 19, 23, 47], zone="region_static")
d.blocks.set(11, [ 16, 36, 37, 20, 52, 38, 39, 53], zone="region_static")
d.blocks.set(12, [ 36, 17, 21, 37, 38, 54, 55, 39], zone="region_dynamic")
d.blocks.set(13, [ 17, 18, 22, 21, 54, 56, 57, 55], zone="region_dynamic")
d.blocks.set(14, [ 18, 46, 47, 22, 56, 48, 49, 57], zone="region_dynamic")
d.blocks.set(15, [ 46, 19, 23, 47, 48, 58, 59, 49], zone="region_static")
d.blocks.set(16, [ 52, 38, 39, 53, 24, 40, 41, 28], zone="region_static")
d.blocks.set(17, [ 38, 54, 55, 39, 40, 25, 29, 41], zone="region_static")
d.blocks.set(18, [ 54, 56, 57, 55, 25, 26, 30, 29], zone="region_static")
d.blocks.set(19, [ 56, 48, 49, 57, 26, 50, 51, 30], zone="region_static")
d.blocks.set(20, [ 48, 58, 59, 49, 50, 27, 31, 51], zone="region_static")

d.blocks.distribution.set( 8, [70, 1, 20])

d.blocks.distribution.set( 6, "x",100)
d.blocks.distribution.set( 7, "x", 10)
d.blocks.distribution.set( 9, "x", 10)
d.blocks.distribution.set(10, "x",100)

d.blocks.distribution.set( 3, "z", 60)
d.blocks.distribution.set(13, "z", 30)
d.blocks.distribution.set(18, "z", 90)

d.boundaryFaces.set(1, "front", [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, \
                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "y-")
d.boundaryFaces.set(2, "back",  [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, \
                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "y+")
d.boundaryFaces.set(3, "infinity", [  1,  2,  3,  4,  5], "z-")
d.boundaryFaces.set(4, "infinity", [ 16, 17, 18, 19, 20], "z+")
d.boundaryFaces.set(5, "infinity", [  1,  6, 11, 16], "x-")
d.boundaryFaces.set(6, "infinity", [  5, 10, 15, 20], "x+")

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

        d.manual("// This is just to demonstrate the manual write")
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
