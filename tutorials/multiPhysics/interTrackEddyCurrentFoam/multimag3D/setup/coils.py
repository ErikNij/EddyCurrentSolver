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

from foamTools.coils import (inductorCoils, writeCoilFeatureEdgeMeshes,
                             writeEdgeBiotSavartProperties, writeFrequency)

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

# --------------------------------------------------------------------------- #
# --- Coils ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coils = inductorCoils(par.coil_setup, "coil", par.coil_bundle, par.coil_path,
                      par.coils_current, par.coils_n, par.coils_step,
                      reverse=par.coils_reverse, origin=par.coils_origin,
                      axis=2, scale=par.coil_scale, period=par.coils_period)

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    writeCoilFeatureEdgeMeshes(par.dir_case, coils)
    writeEdgeBiotSavartProperties(par.dir_case, coils, par.coils_nNonOrto)
    writeFrequency(par.dir_case, par.coils_frequency)

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
