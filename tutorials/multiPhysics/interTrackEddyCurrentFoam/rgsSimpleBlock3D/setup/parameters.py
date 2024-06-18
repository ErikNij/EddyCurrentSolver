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

geo_scale       = 1e-3

geo_angle       = 0.0
geo_alpha       = geo_angle/180.0 * m.pi

geo_X           = dict()
geo_X["C"]      = 35.0

geo_Y           = dict()
geo_Y["C"]      = 75.0
geo_Y["C"]      = 75.0

geo_Z           = dict()
geo_Z["C"]      = 20.0

geo_Y["shift"]  = geo_Z["C"]/2.0 * m.tan(geo_alpha)

# --------------------------------------------------------------------------- #
# --- Mesh ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

mesh_scale      = 1.0
mesh_space      = 4.0

mesh_normal     = -1

mesh_X          = dict()
mesh_X["inf"]   = mesh_space * 1.5 * geo_X["C"]

mesh_Y          = dict()
mesh_Y["inf"]   = mesh_space * geo_Y["C"]

mesh_Z          = dict()
mesh_Z["inf+"]  = mesh_space * 3.0 * geo_Z["C"]/2.0 + geo_Z["C"]
mesh_Z["inf-"]  = mesh_space * 3.0 * geo_Z["C"]/2.0

mesh            = {"normal": mesh_normal}

# --------------------------------------------------------------------------- #
# --- Coils ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coil_scale      = 1e-3

# --------------------------------------------------------------------------- #

#coil_bundle     = {"shape": "rectangle",
coil_bundle     = {"shape": "point",
                    "n":     4,
                    "r":     5.0,
                    "z":     5.0}

#coil_path       = {"shape": "rectangle",
coil_path       = {"shape": "racetrack",
                    "n":     12,
                    "r":     8.0,
                    "x":     65.0,
                    "y":     120.0}

coils_n         = 3
coils_step      = 13.0
coils_origin    = [0.0, 0.0, 9.0]

coils_current   = 1000.0 * m.sqrt(2)
coils_nNonOrto  = 10
coils_frequency = 1000.0

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
