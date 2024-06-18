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

geo_scale      = 1e-3

geo_R          = dict()
geo_R[0]       =  53.0
geo_R[1]       =  62.0
geo_R[2]       =  70.0

geo_Z          = dict()
geo_Z[0]       = -34.0
geo_Z[1]       =  -8.0
geo_Z[2]       =   0.0
geo_Z[3]       =  22.0
geo_Z[4]       =  27.095
geo_Z[5]       =  30.0
geo_Z[6]       =  70.0
geo_Z[7]       = 142.0

geo_alpha      = m.atan(geo_Z[5]/geo_R[0])                # Cone angle (rad)

geo_magG       = (geo_Z[2] - geo_Z[1]) * m.cos(geo_alpha) # Non-conducting gap size
geo_G          = geo_magG * np.array([m.sin(geo_alpha), -m.cos(geo_alpha)])

geo_Z["solid"] =  geo_Z[3]                                # Solid height
geo_R["solid"] =  geo_Z[3] / m.tan(geo_alpha)             # Solid wall contact radius

# --------------------------------------------------------------------------- #
# --- Mesh ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

mesh_scale     = 1.0
mesh_space     = 4.0

mesh_normal    = 1
mesh_thickness = 10.0
mesh_angle     =  5.0

mesh_R         = dict()
mesh_R["axis"] = 0.01                               # Axis patch radius
mesh_R["inf"]  = mesh_space * geo_R[2]              # Infinity patch radius

mesh_Z         = dict()
mesh_Z["C"]    = geo_Z[0] + 0.5*(geo_Z[7] - geo_Z[0])
mesh_Z["inf+"] = mesh_space * abs(geo_Z[7] - mesh_Z["C"]) # Infinity patch radius +
mesh_Z["inf-"] = mesh_space * abs(geo_Z[0] - mesh_Z["C"]) # Infinity patch radius -
mesh_Z["inf"]  = max(mesh_Z["inf+"], mesh_Z["inf-"])      # Infinity patch radius

mesh_R["inf"]  = max(mesh_R["inf"], mesh_Z["inf"])
mesh_Z["inf"]  = max(mesh_Z["inf"], mesh_R["inf"])

mesh           = {"normal": mesh_normal,
                  "wedge": False,
                  "extent": mesh_thickness}

# --------------------------------------------------------------------------- #
# --- Coils ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coil_scale      = 1e-3

coil_bundle     = {"shape": "rectangle",
                   "n":     10,
                   "r":     10.0,
                   "z":     8.0}

coil_path       = {"shape": "loop",
                   "n":     36,
                   "r":     100.0 + coil_bundle["r"]/2.0}

coils_n         = 10
coils_step      = 14.9
coils_origin    = [0.0, 0.0, 3.0 + coil_bundle["z"]/2.0]

coils_current   = m.sqrt(2.0) * 260.0
coils_nNonOrto  = 10
coils_frequency = 6300.0

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
