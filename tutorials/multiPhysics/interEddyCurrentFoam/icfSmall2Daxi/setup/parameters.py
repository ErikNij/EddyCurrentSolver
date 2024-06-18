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
geo_R[0]       =  85.0
geo_R[1]       = geo_R[0] + 1.0
geo_R[2]       = geo_R[1] + 14.0
geo_R[3]       = geo_R[2] + 15.0

geo_Z          = dict()
geo_Z[0]       =  -1.0
geo_Z[1]       =   0.0
geo_Z[2]       = 235.0
geo_Z[3]       = 383.0

# --------------------------------------------------------------------------- #
# --- Mesh ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

mesh_scale     = 2.0
mesh_space     = 4.0

mesh_normal    = 1
mesh_thickness = 10.0
mesh_angle     =  5.0

mesh_R         = dict()
mesh_R["axis"] = 0.0                                # Axis patch radius
mesh_R["inf"]  = mesh_space * geo_R[1]              # Infinity patch radius
#mesh_R["inf"]  = mesh_space * geo_R[3]              # Infinity patch radius (j0)

mesh_Z         = dict()
mesh_Z["C"]    = geo_Z[0] + 0.5*(geo_Z[3] - geo_Z[0])
mesh_Z["inf+"] = mesh_space * abs(geo_Z[3] - mesh_Z["C"]) # Infinity patch radius +
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

coil_bundle     = {"shape": "point",
#coil_bundle     = {"shape": "rectangle",
                   "n":     10,
                   "r":     (geo_R[3] - geo_R[2]),
                   "z":     15.0}

coil_path       = {"shape": "loop",
                   "n":     36,
                   #"r":     geo_R[2] + coil_bundle["r"]/2.0}
                   "r":     geo_R[2]}

coils_n         = 12
coils_step      = 5.0 + coil_bundle["z"]
#coils_step      = 5.0 + coil_bundle["z"] -1.363636
coils_origin    = [0.0, 0.0, geo_Z[1] + coil_bundle["z"]/2.0]
#coils_origin    = [0.0, 0.0, geo_Z[1] + coil_bundle["z"]]

#coils_current   = m.sqrt(2.0) * 1750.0 (fL ~ 2.9e+5)
#coils_current   = m.sqrt(2.0) * 2020.0
#coils_current   = m.sqrt(2.0) * 2260.0
coils_current   = m.sqrt(2.0) * 2460.0
coils_nNonOrto  = 10
coils_frequency = 330.0

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
