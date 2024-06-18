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

# Basic

geo_scale      =   1.0e-3

geo_rad        = m.pi / 180.0
geo_cos45      = m.cos(45*geo_rad)
geo_sin45      = m.sin(45*geo_rad)

geo_f001       =     0.25
geo_f112       =   1.0

geo_melt_radius_f   =   1.01
geo_coil_radius_f   =   1.05

# Melt

geo_melt_r     =   6.0
geo_melt_x     =  35.0
geo_melt_y     =  75.0

# Frame

geo_frame_d    =  10.0
geo_frame_x    = geo_melt_x + geo_frame_d
geo_frame_y    = geo_melt_y + geo_frame_d

# Substrate

geo_subs_x     =  78.0
geo_subs_y     =  79.0
geo_subs_wafer = geo_subs_y - geo_melt_y

# Coil
geo_coil_x     =  65.0
geo_coil_y     = 120.0
geo_coil_d     =   5.0
geo_coil_r     =   8.0

# Geometry x

geo_x00        = geo_melt_x - geo_melt_r
geo_x01        = geo_melt_x - geo_melt_r/2.0
geo_x02        = geo_melt_x
geo_x03        = geo_melt_x + geo_subs_wafer
geo_x04        = geo_melt_x + geo_frame_d
geo_x05        = geo_coil_x - geo_coil_r
geo_x06        = geo_coil_x - geo_coil_d/2.0
geo_x07        = geo_coil_x + geo_coil_d/2.0
geo_x08        = geo_coil_x + geo_coil_r
geo_x09        = geo_subs_x
geo_x10        = 150.0

geo_xp09       = geo_x00 + (geo_f001 + (1.0-geo_f001)*geo_cos45)*geo_melt_r/2.0
geo_xp13       = geo_x00 + geo_cos45*geo_melt_r*geo_melt_radius_f
geo_xp16       = geo_x01
geo_xp18       = geo_x00 + (geo_f112 + (1.0-geo_f112)*geo_cos45)*(geo_melt_r+geo_subs_wafer)
geo_xp67       = geo_x05 + geo_cos45*(geo_coil_r-geo_coil_d/2.0)*geo_coil_radius_f
geo_xp73       = geo_x05 + geo_cos45*(geo_coil_r+geo_coil_d/2.0)*geo_coil_radius_f
geo_xp79       = geo_x05 + geo_cos45*2.0*geo_coil_r*geo_coil_radius_f

# Geometry y

geo_y00        = geo_melt_y - geo_melt_r
geo_y01        = geo_melt_y - geo_melt_r/2.0
geo_y02        = geo_melt_y
geo_y03        = geo_melt_y + geo_subs_wafer
geo_y04        = geo_melt_y + geo_frame_d
geo_y05        = geo_coil_y - geo_coil_r
geo_y06        = geo_coil_y - geo_coil_d/2.0
geo_y07        = geo_coil_y + geo_coil_d/2.0
geo_y08        = geo_coil_y + geo_coil_r
geo_y09        = geo_coil_y + (geo_subs_x - geo_coil_x)
geo_y10        = 200.0

geo_yp09       = geo_y00 + (geo_f001 + (1.0-geo_f001)*geo_sin45)*geo_melt_r/2.0
geo_yp13       = geo_y00 + geo_sin45*geo_melt_r*geo_melt_radius_f
geo_yp18       = geo_y00 + (geo_f112 + (1.0-geo_f112)*geo_sin45)*(geo_melt_r+geo_subs_wafer)
geo_yp19       = geo_melt_y - geo_subs_wafer
geo_yp67       = geo_y05 + geo_sin45*(geo_coil_r-geo_coil_d/2.0)*geo_coil_radius_f
geo_yp73       = geo_y05 + geo_sin45*(geo_coil_r+geo_coil_d/2.0)*geo_coil_radius_f
geo_yp79       = geo_y05 + geo_sin45*2.0*geo_coil_r*geo_coil_radius_f

# Geometry z

geo_z00        = -60.0
geo_z01        =  -5.0
geo_z02        =  -2.0 # NOTE: The block between 0.0 and -2-0 will be split in half to get level at -1.0
geo_z03        =   0.0
geo_z04        =   7.0
geo_z05        =  12.0
geo_z06        =  20.0
geo_z07        =  25.0
geo_z08        =  33.0
geo_z09        =  38.0
geo_z10        =  50.0
geo_z11        =  90.0
geo_z12        = 140.0

# Grading and Distribution

geo_grade_ref = 2.0

geo_dist_ref = geo_frame_d/2.0

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

# Vertices

d.vertices.set(  0, [      0.0,      0.0, geo_z00])
d.vertices.set(  1, [  geo_x00,      0.0, geo_z00])
d.vertices.set(  2, [  geo_x00,  geo_y00, geo_z00])
d.vertices.set(  3, [      0.0,  geo_y00, geo_z00])
# WARNING: Vertice 4 and 5 are missing. Remember this!
d.vertices.set(  6, [  geo_x00,  geo_y01, geo_z00])
d.vertices.set(  7, [      0.0,  geo_y01, geo_z00])
d.vertices.set(  8, [  geo_x01,  geo_y00, geo_z00])
d.vertices.set(  9, [ geo_xp09, geo_yp09, geo_z00])
d.vertices.set( 10, [  geo_x01,      0.0, geo_z00])
d.vertices.set( 11, [  geo_x00,  geo_y02, geo_z00])
d.vertices.set( 12, [      0.0,  geo_y02, geo_z00])
d.vertices.set( 13, [ geo_xp13, geo_yp13, geo_z00])
d.vertices.set( 14, [  geo_x02,  geo_y00, geo_z00])
d.vertices.set( 15, [  geo_x02,      0.0, geo_z00])
d.vertices.set( 16, [ geo_xp16,  geo_y03, geo_z00])
d.vertices.set( 17, [      0.0,  geo_y03, geo_z00])
d.vertices.set( 18, [ geo_xp18, geo_yp18, geo_z00])
d.vertices.set( 19, [  geo_x03, geo_yp19, geo_z00])
d.vertices.set( 20, [  geo_x03,      0.0, geo_z00])
d.vertices.set( 21, [ geo_xp16,  geo_y04, geo_z00])
d.vertices.set( 22, [      0.0,  geo_y04, geo_z00])
d.vertices.set( 23, [ geo_xp18,  geo_y04, geo_z00])
d.vertices.set( 24, [  geo_x04, geo_yp18, geo_z00])
d.vertices.set( 25, [  geo_x04,  geo_y04, geo_z00])
d.vertices.set( 26, [  geo_x04, geo_yp19, geo_z00])
d.vertices.set( 27, [  geo_x04,      0.0, geo_z00])
d.vertices.set( 28, [ geo_xp16,  geo_y05, geo_z00])
d.vertices.set( 29, [      0.0,  geo_y05, geo_z00])
d.vertices.set( 30, [ geo_xp18,  geo_y05, geo_z00])
d.vertices.set( 31, [  geo_x04,  geo_y05, geo_z00])

d.vertices.set( 32, [  geo_x05,  geo_y04, geo_z00])
d.vertices.set( 33, [  geo_x05,  geo_y05, geo_z00])
d.vertices.set( 34, [  geo_x05, geo_yp18, geo_z00])
d.vertices.set( 35, [  geo_x05, geo_yp19, geo_z00])
d.vertices.set( 36, [  geo_x05,      0.0, geo_z00])

d.vertices.set( 37, [  geo_x06,  geo_y05, geo_z00])
d.vertices.set( 38, [  geo_x06,  geo_y04, geo_z00])
d.vertices.set( 39, [  geo_x06, geo_yp18, geo_z00])
d.vertices.set( 40, [  geo_x06, geo_yp19, geo_z00])
d.vertices.set( 41, [  geo_x06,      0.0, geo_z00])

d.vertices.set( 42, [  geo_x07,  geo_y05, geo_z00])
d.vertices.set( 43, [  geo_x07,  geo_y04, geo_z00])
d.vertices.set( 44, [  geo_x07, geo_yp18, geo_z00])
d.vertices.set( 45, [  geo_x07, geo_yp19, geo_z00])
d.vertices.set( 46, [  geo_x07,      0.0, geo_z00])

d.vertices.set( 47, [  geo_x08,  geo_y05, geo_z00])
d.vertices.set( 48, [  geo_x08,  geo_y04, geo_z00])
d.vertices.set( 49, [  geo_x08, geo_yp18, geo_z00])
d.vertices.set( 50, [  geo_x08, geo_yp19, geo_z00])
d.vertices.set( 51, [  geo_x08,      0.0, geo_z00])

d.vertices.set( 52, [  geo_x09,  geo_y05, geo_z00])
d.vertices.set( 53, [  geo_x09,  geo_y04, geo_z00])
d.vertices.set( 54, [  geo_x09, geo_yp18, geo_z00])
d.vertices.set( 55, [  geo_x09, geo_yp19, geo_z00])
d.vertices.set( 56, [  geo_x09,      0.0, geo_z00])

d.vertices.set( 57, [  geo_x10,  geo_y05, geo_z00])
d.vertices.set( 58, [  geo_x10,  geo_y04, geo_z00])
d.vertices.set( 59, [  geo_x10, geo_yp18, geo_z00])
d.vertices.set( 60, [  geo_x10, geo_yp19, geo_z00])
d.vertices.set( 61, [  geo_x10,      0.0, geo_z00])

d.vertices.set( 62, [      0.0,  geo_y06, geo_z00])
d.vertices.set( 63, [ geo_xp16,  geo_y06, geo_z00])
d.vertices.set( 64, [ geo_xp18,  geo_y06, geo_z00])
d.vertices.set( 65, [  geo_x04,  geo_y06, geo_z00])
d.vertices.set( 66, [  geo_x05,  geo_y06, geo_z00])
d.vertices.set( 67, [ geo_xp67, geo_yp67, geo_z00])

d.vertices.set( 68, [      0.0,  geo_y07, geo_z00])
d.vertices.set( 69, [ geo_xp16,  geo_y07, geo_z00])
d.vertices.set( 70, [ geo_xp18,  geo_y07, geo_z00])
d.vertices.set( 71, [  geo_x04,  geo_y07, geo_z00])
d.vertices.set( 72, [  geo_x05,  geo_y07, geo_z00])
d.vertices.set( 73, [ geo_xp73, geo_yp73, geo_z00])

d.vertices.set( 74, [      0.0,  geo_y08, geo_z00])
d.vertices.set( 75, [ geo_xp16,  geo_y08, geo_z00])
d.vertices.set( 76, [ geo_xp18,  geo_y08, geo_z00])
d.vertices.set( 77, [  geo_x04,  geo_y08, geo_z00])
d.vertices.set( 78, [  geo_x05,  geo_y08, geo_z00])
d.vertices.set( 79, [ geo_xp79, geo_yp79, geo_z00])
d.vertices.set( 80, [  geo_x09, geo_yp79, geo_z00])
d.vertices.set( 81, [  geo_x10, geo_yp79, geo_z00])

d.vertices.set( 82, [      0.0,  geo_y09, geo_z00])
d.vertices.set( 83, [ geo_xp16,  geo_y09, geo_z00])
d.vertices.set( 84, [ geo_xp18,  geo_y09, geo_z00])
d.vertices.set( 85, [  geo_x04,  geo_y09, geo_z00])
d.vertices.set( 86, [  geo_x05,  geo_y09, geo_z00])
d.vertices.set( 87, [ geo_xp79,  geo_y09, geo_z00])
d.vertices.set( 88, [  geo_x09,  geo_y09, geo_z00])
d.vertices.set( 89, [  geo_x10,  geo_y09, geo_z00])

d.vertices.set( 90, [      0.0,  geo_y10, geo_z00])
d.vertices.set( 91, [ geo_xp16,  geo_y10, geo_z00])
d.vertices.set( 92, [ geo_xp18,  geo_y10, geo_z00])
d.vertices.set( 93, [  geo_x04,  geo_y10, geo_z00])
d.vertices.set( 94, [  geo_x05,  geo_y10, geo_z00])
d.vertices.set( 95, [ geo_xp79,  geo_y10, geo_z00])
d.vertices.set( 96, [  geo_x09,  geo_y10, geo_z00])
d.vertices.set( 97, [  geo_x10,  geo_y10, geo_z00])

baseVertices = [ i for i in range(4) ] + [ i+6 for i in range(len(d.vertices.labels)-4) ]

d.vertices.copyTranslate( 100, baseVertices, [0.0, 0.0, geo_z01-geo_z00])
d.vertices.copyTranslate( 200, baseVertices, [0.0, 0.0, geo_z02-geo_z00])
d.vertices.copyTranslate( 300, baseVertices, [0.0, 0.0, geo_z03-geo_z00])
d.vertices.copyTranslate( 400, baseVertices, [0.0, 0.0, geo_z04-geo_z00])
d.vertices.copyTranslate( 500, baseVertices, [0.0, 0.0, geo_z05-geo_z00])
d.vertices.copyTranslate( 600, baseVertices, [0.0, 0.0, geo_z06-geo_z00])
d.vertices.copyTranslate( 700, baseVertices, [0.0, 0.0, geo_z07-geo_z00])
d.vertices.copyTranslate( 800, baseVertices, [0.0, 0.0, geo_z08-geo_z00])
d.vertices.copyTranslate( 900, baseVertices, [0.0, 0.0, geo_z09-geo_z00])
d.vertices.copyTranslate(1000, baseVertices, [0.0, 0.0, geo_z10-geo_z00])
d.vertices.copyTranslate(1100, baseVertices, [0.0, 0.0, geo_z11-geo_z00])
d.vertices.copyTranslate(1200, baseVertices, [0.0, 0.0, geo_z12-geo_z00])

# Blocks

d.blocks.set(  0, [  0,   1,   2,   3, 100, 101, 102, 103], zone="fluid")
d.blocks.set(  1, [  3,   2,   6,   7, 103, 102, 106, 107], zone="fluid")
d.blocks.set(  2, [  2,   8,   9,   6, 102, 108, 109, 106], zone="fluid")
d.blocks.set(  3, [  1,  10,   8,   2, 101, 110, 108, 102], zone="fluid")
d.blocks.set(  4, [  7,   6,  11,  12, 107, 106, 111, 112], zone="fluid")
d.blocks.set(  5, [  6,   9,  13,  11, 106, 109, 113, 111], zone="fluid")
d.blocks.set(  6, [  8,  14,  13,   9, 108, 114, 113, 109], zone="fluid")
d.blocks.set(  7, [ 10,  15,  14,   8, 110, 115, 114, 108], zone="fluid")
d.blocks.set(  8, [ 12,  11,  16,  17, 112, 111, 116, 117], zone="frame")
d.blocks.set(  9, [ 11,  13,  18,  16, 111, 113, 118, 116], zone="frame")
d.blocks.set( 10, [ 14,  19,  18,  13, 114, 119, 118, 113], zone="frame")
d.blocks.set( 11, [ 15,  20,  19,  14, 115, 120, 119, 114], zone="frame")
d.blocks.set( 12, [ 17,  16,  21,  22, 117, 116, 121, 122], zone="frame")
d.blocks.set( 13, [ 16,  18,  23,  21, 116, 118, 123, 121], zone="frame")
d.blocks.set( 14, [ 18,  24,  25,  23, 118, 124, 125, 123], zone="frame")
d.blocks.set( 15, [ 19,  26,  24,  18, 119, 126, 124, 118], zone="frame")
d.blocks.set( 16, [ 20,  27,  26,  19, 120, 127, 126, 119], zone="frame")
d.blocks.set( 17, [ 22,  21,  28,  29, 122, 121, 128, 129], zone="background")
d.blocks.set( 18, [ 21,  23,  30,  28, 121, 123, 130, 128], zone="background")
d.blocks.set( 19, [ 23,  25,  31,  30, 123, 125, 131, 130], zone="background")

d.blocks.set( 20, [ 25,  32,  33,  31, 125, 132, 133, 131], zone="background")
d.blocks.set( 21, [ 24,  34,  32,  25, 124, 134, 132, 125], zone="background")
d.blocks.set( 22, [ 26,  35,  34,  24, 126, 135, 134, 124], zone="background")
d.blocks.set( 23, [ 27,  36,  35,  26, 127, 136, 135, 126], zone="background")

d.blocks.set( 24, [ 32,  38,  37,  33, 132, 138, 137, 133], zone="background")
d.blocks.set( 25, [ 34,  39,  38,  32, 134, 139, 138, 132], zone="background")
d.blocks.set( 26, [ 35,  40,  39,  34, 135, 140, 139, 134], zone="background")
d.blocks.set( 27, [ 36,  41,  40,  35, 136, 141, 140, 135], zone="background")

d.blocks.set( 28, [ 38,  43,  42,  37, 138, 143, 142, 137], zone="coil")
d.blocks.set( 29, [ 39,  44,  43,  38, 139, 144, 143, 138], zone="coil")
d.blocks.set( 30, [ 40,  45,  44,  39, 140, 145, 144, 139], zone="coil")
d.blocks.set( 31, [ 41,  46,  45,  40, 141, 146, 145, 140], zone="coil")

d.blocks.set( 32, [ 43,  48,  47,  42, 143, 148, 147, 142], zone="background")
d.blocks.set( 33, [ 44,  49,  48,  43, 144, 149, 148, 143], zone="background")
d.blocks.set( 34, [ 45,  50,  49,  44, 145, 150, 149, 144], zone="background")
d.blocks.set( 35, [ 46,  51,  50,  45, 146, 151, 150, 145], zone="background")

d.blocks.set( 36, [ 48,  53,  52,  47, 148, 153, 152, 147], zone="background")
d.blocks.set( 37, [ 49,  54,  53,  48, 149, 154, 153, 148], zone="background")
d.blocks.set( 38, [ 50,  55,  54,  49, 150, 155, 154, 149], zone="background")
d.blocks.set( 39, [ 51,  56,  55,  50, 151, 156, 155, 150], zone="background")

d.blocks.set( 40, [ 53,  58,  57,  52, 153, 158, 157, 152], zone="background")
d.blocks.set( 41, [ 54,  59,  58,  53, 154, 159, 158, 153], zone="background")
d.blocks.set( 42, [ 55,  60,  59,  54, 155, 160, 159, 154], zone="background")
d.blocks.set( 43, [ 56,  61,  60,  55, 156, 161, 160, 155], zone="background")

d.blocks.set( 44, [ 29,  28,  63,  62, 129, 128, 163, 162], zone="background")
d.blocks.set( 45, [ 28,  30,  64,  63, 128, 130, 164, 163], zone="background")
d.blocks.set( 46, [ 30,  31,  65,  64, 130, 131, 165, 164], zone="background")
d.blocks.set( 47, [ 31,  33,  66,  65, 131, 133, 166, 165], zone="background")
d.blocks.set( 48, [ 33,  37,  67,  66, 133, 137, 167, 166], zone="background")

d.blocks.set( 49, [ 62,  63,  69,  68, 162, 163, 169, 168], zone="coil")
d.blocks.set( 50, [ 63,  64,  70,  69, 163, 164, 170, 169], zone="coil")
d.blocks.set( 51, [ 64,  65,  71,  70, 164, 165, 171, 170], zone="coil")
d.blocks.set( 52, [ 65,  66,  72,  71, 165, 166, 172, 171], zone="coil")
d.blocks.set( 53, [ 66,  67,  73,  72, 166, 167, 173, 172], zone="coil")
d.blocks.set( 54, [ 37,  42,  73,  67, 137, 142, 173, 167], zone="coil")

d.blocks.set( 55, [ 68,  69,  75,  74, 168, 169, 175, 174], zone="background")
d.blocks.set( 56, [ 69,  70,  76,  75, 169, 170, 176, 175], zone="background")
d.blocks.set( 57, [ 70,  71,  77,  76, 170, 171, 177, 176], zone="background")
d.blocks.set( 58, [ 71,  72,  78,  77, 171, 172, 178, 177], zone="background")
d.blocks.set( 59, [ 72,  73,  79,  78, 172, 173, 179, 178], zone="background")
d.blocks.set( 60, [ 42,  47,  79,  73, 142, 147, 179, 173], zone="background")
d.blocks.set( 61, [ 47,  52,  80,  79, 147, 152, 180, 179], zone="background")
d.blocks.set( 62, [ 52,  57,  81,  80, 152, 157, 181, 180], zone="background")

d.blocks.set( 63, [ 74,  75,  83,  82, 174, 175, 183, 182], zone="background")
d.blocks.set( 64, [ 75,  76,  84,  83, 175, 176, 184, 183], zone="background")
d.blocks.set( 65, [ 76,  77,  85,  84, 176, 177, 185, 184], zone="background")
d.blocks.set( 66, [ 77,  78,  86,  85, 177, 178, 186, 185], zone="background")
d.blocks.set( 67, [ 78,  79,  87,  86, 178, 179, 187, 186], zone="background")
d.blocks.set( 68, [ 79,  80,  88,  87, 179, 180, 188, 187], zone="background")
d.blocks.set( 69, [ 80,  81,  89,  88, 180, 181, 189, 188], zone="background")

d.blocks.set( 70, [ 82,  83,  91,  90, 182, 183, 191, 190], zone="background")
d.blocks.set( 71, [ 83,  84,  92,  91, 183, 184, 192, 191], zone="background")
d.blocks.set( 72, [ 84,  85,  93,  92, 184, 185, 193, 192], zone="background")
d.blocks.set( 73, [ 85,  86,  94,  93, 185, 186, 194, 193], zone="background")
d.blocks.set( 74, [ 86,  87,  95,  94, 186, 187, 195, 194], zone="background")
d.blocks.set( 75, [ 87,  88,  96,  95, 187, 188, 196, 195], zone="background")
d.blocks.set( 76, [ 88,  89,  97,  96, 188, 189, 197, 196], zone="background")

baseBlocks = [ i for i in range(len(d.blocks.labels)) ]

d.blocks.copyShiftVerticeLabels(  100, baseBlocks,  100)
d.blocks.copyShiftVerticeLabels(  200, baseBlocks,  200)
d.blocks.copyShiftVerticeLabels(  300, baseBlocks,  300)
d.blocks.copyShiftVerticeLabels(  400, baseBlocks,  400)
d.blocks.copyShiftVerticeLabels(  500, baseBlocks,  500)
d.blocks.copyShiftVerticeLabels(  600, baseBlocks,  600)
d.blocks.copyShiftVerticeLabels(  700, baseBlocks,  700)
d.blocks.copyShiftVerticeLabels(  800, baseBlocks,  800)
d.blocks.copyShiftVerticeLabels(  900, baseBlocks,  900)
d.blocks.copyShiftVerticeLabels( 1000, baseBlocks, 1000)
d.blocks.copyShiftVerticeLabels( 1100, baseBlocks, 1100)

# Gradings

d.blocks.grading.set(   76, [geo_grade_ref, geo_grade_ref, 1.0/geo_grade_ref])
d.blocks.grading.set( 1176, [geo_grade_ref, geo_grade_ref,     geo_grade_ref])

# Distributions

d.blocks.distribution.set(    0, "x", max(1, int(m.floor(geo_melt_x/geo_dist_ref))))
d.blocks.distribution.set(   20, "x", max(1, int(m.floor((geo_x05-geo_x04)/geo_dist_ref))))
d.blocks.distribution.set(   76, "x", max(1, int(m.floor((geo_x10-geo_x09)/geo_dist_ref/geo_grade_ref))))

d.blocks.distribution.set(    0, "y", max(1, int(m.floor(geo_melt_y/geo_dist_ref))))
d.blocks.distribution.set(   20, "y", max(1, int(m.floor((geo_y05-geo_y04)/geo_dist_ref))))
d.blocks.distribution.set(   76, "y", max(1, int(m.floor((geo_y10-geo_y09)/geo_dist_ref/geo_grade_ref))))

d.blocks.distribution.set(    0, "z", max(1, int(m.floor((geo_z01-geo_z00)/geo_dist_ref/geo_grade_ref))))
d.blocks.distribution.set(  100, "z", max(1, int(m.floor((geo_z02-geo_z01)/geo_dist_ref))))
d.blocks.distribution.set(  200, "z", max(1, int(m.floor((geo_z03-geo_z02)/geo_dist_ref))))
d.blocks.distribution.set(  300, "z", max(1, int(m.floor((geo_z04-geo_z03)/geo_dist_ref))))
d.blocks.distribution.set(  400, "z", max(1, int(m.floor((geo_z05-geo_z04)/geo_dist_ref))))
d.blocks.distribution.set(  500, "z", max(1, int(m.floor((geo_z06-geo_z05)/geo_dist_ref))))
d.blocks.distribution.set(  600, "z", max(1, int(m.floor((geo_z07-geo_z06)/geo_dist_ref))))
d.blocks.distribution.set(  700, "z", max(1, int(m.floor((geo_z08-geo_z07)/geo_dist_ref))))
d.blocks.distribution.set(  800, "z", max(1, int(m.floor((geo_z09-geo_z08)/geo_dist_ref))))
d.blocks.distribution.set(  900, "z", max(1, int(m.floor((geo_z10-geo_z09)/geo_dist_ref))))
d.blocks.distribution.set( 1000, "z", max(1, int(m.floor((geo_z11-geo_z10)/geo_dist_ref))))
d.blocks.distribution.set( 1100, "z", max(1, int(m.floor((geo_z12-geo_z11)/geo_dist_ref/geo_grade_ref))))

## Boundary faces

#d.boundaryFaces.set(0, "mirror_x", [0, 1, 4, 8, 12, 17], "x-")
#d.boundaryFaces.set(3, "mirror_y", [0, 3, 7, 11, 16, 23], "y-")
## ...

#d.boundaryFaces.set(1, "infinity", [24, 25, 26, 27], "x+")
#d.boundaryFaces.set(4, "infinity", [17, 18, 19, 20], "y+")
## ...

#d.boundaryFaces.set(5, "infinity", [ i for i in range(24) ], "z-")
#d.boundaryFaces.set(6, "infinity", [ i+600 for i in range(24) ], "z+")

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

    if d.boundarySubDict("mirror_x", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_y", "patch"):

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
