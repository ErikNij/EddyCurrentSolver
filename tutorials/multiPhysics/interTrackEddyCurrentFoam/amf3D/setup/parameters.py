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

# --------------------------------------------------------------------------- #
# --- Geometry -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

geo_scale  = 1e-3
geo_space  = 4.0

R = 89.5
H = 59.0

geo_r0     =   R * 0.85
geo_r1     =   R
geo_r2     =   R * geo_space

geo_z0     =  -H * geo_space
geo_z1     =  -H
geo_z2     =  -H + 0.15 * R
geo_z3     =   H - 0.15 * R
geo_z4     =   H
geo_z5     =   H * geo_space

# --------------------------------------------------------------------------- #
# --- Mesh ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

mesh_scale =    0.25

mesh_normal    = -1
mesh_f     =    1.0

mesh           = {"normal": mesh_normal}

# --------------------------------------------------------------------------- #
# --- Coils ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coil_scale      = 1e-3

#coil_bundle     = {"shape": "rectangle",
coil_bundle     = {"shape": "point",
                    "n":     5,
                    "r":     8.0,
                    "z":     8.0}

coil_path       = {"shape": "racetrack",
                    "n":     9,
                    "r":     20.0 + coil_bundle["r"]/2.0,
                    "x":     164.0 + coil_bundle["r"]/2.0,
                    "y":     164.0 + coil_bundle["r"]/2.0}

coils_n         = 24
coils_step      = 9.25
coils_origin    = [0.0, 0.0, - 106.375]

coils_current   = 10.345
coils_reverse   = False
coils_nNonOrto  = 10
coils_frequency = 96.11
#coils_frequency = 480.58

# --------------------------------------------------------------------------- #
# --- Directories ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

dir_case = os.path.realpath(__dir__ + "/" + "..")

dir_setup = os.path.realpath(dir_case + "/" + "setup")
dir_results = os.path.realpath(dir_case + "/" + "results")

dir_0 = os.path.realpath(dir_case + "/" + "0")
dir_constant = os.path.realpath(dir_case + "/" + "constant")
dir_system = os.path.realpath(dir_case + "/" + "system")

dir_polyMesh = os.path.realpath(dir_constant + "/" + "polyMesh")
dir_featureEdgeMesh = os.path.realpath(dir_constant + "/" + "featureEdgeMesh")
dir_triSurface = os.path.realpath(dir_constant + "/" + "triSurface")

for d in [dir_polyMesh, dir_triSurface, dir_featureEdgeMesh]:
    if not os.path.exists(d): os.makedirs(d)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
