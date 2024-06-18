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

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lscm
from matplotlib.cm import _reverse_cmap_spec as rcms
from cycler import cycler as ccycler

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def colors():
    """Set up HZDR colors.
    """

    def h2c(h): return mpl.colors.hex2color(h)
    def c2h(c): return mpl.colors.rgb2hex(c)

    c = mpl.colors.ColorConverter.colors
    d = mpl.cm.cmap_d

    # Primary colors
    c["hzdr-blue"]      = h2c("#00589C")
    c["hzdr-orange"]    = h2c("#CF6800")
    c["hzdr-darkblue"]  = h2c("#003E6E")
    c["hzdr-gray1"]     = h2c("#515151")
    c["hzdr-gray2"]     = h2c("#9C9C9C")
    c["hzdr-gray3"]     = h2c("#B9B9B9")

    # Secondary colors
    c["hzdr-health"]    = h2c("#D42D12") # red
    c["hzdr-struct"]    = h2c("#CF6800") # orange
    c["hzdr-energy"]    = h2c("#E6AF11") # yellow
    c["hzdr-earth"]     = h2c("#144D28") # dark green
    c["hzdr-keytec"]    = h2c("#A9B509") # green
    c["hzdr-aero"]      = h2c("#00A2E0") # cyan

    # Tertiary colors (custom)
    c["hzdr-purple"]    = h2c("#640A78")

    # Alias
    c["hzdr-red"]       = c["hzdr-health"]
    c["hzdr-yellow"]    = c["hzdr-energy"]
    c["hzdr-darkgreen"] = c["hzdr-earth"]
    c["hzdr-green"]     = c["hzdr-keytec"]
    c["hzdr-cyan"]      = c["hzdr-aero"]

    # Create new color cycle
    oldPropCycle = mpl.rcParams["axes.prop_cycle"]
    oldPropStates = dict()
    oldProps = list(oldPropCycle)[0].keys()
    for prop in oldProps:
        oldPropStates[prop] = [cyc[prop] for cyc in list(oldPropCycle)]

    newPropCycle = ccycler("color", ["hzdr-blue", "hzdr-orange", "hzdr-green",
                                     "hzdr-red", "hzdr-yellow",  "hzdr-purple"])

    for prop in oldProps:
        if (prop != "color"):
            newPropCycle += ccycler(prop, oldPropStates[prop])

    # Update default color cycle
    mpl.rcParams["axes.prop_cycle"] = newPropCycle

    # Two-level color maps
    for simple in ["blue", "orange", "green", "red", "yellow", "purple"]:
        cn = "hzdr-"+simple
        d[cn] = lscm.from_list(cn, [(0.0, "white"), (1.0, cn)])
        d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    # Two-level color maps
    cn, cmax = "hzdr-blue-orange", 1.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-blue"),
                                (1/cmax, "hzdr-orange")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    # Three-level color maps
    cn, cmax = "hzdr-blue-white-orange", 2.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-blue"),
                                (1/cmax, "white"),
                                (2/cmax, "hzdr-orange")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-blue-grey3-orange", 2.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-blue"),
                                (1/cmax, "hzdr-gray3"),
                                (2/cmax, "hzdr-orange")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-blue-yellow-red",2.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-blue"),
                                (1/cmax, "hzdr-yellow"),
                                (2/cmax, "hzdr-red")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-orange-green-purple",2.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-orange"),
                                (1/cmax, "hzdr-green"),
                                (2/cmax, "hzdr-purple")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

     # Six-level color maps
    cn, cmax = "hzdr-rainbow1", 5.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-blue"),
                                (1/cmax, "hzdr-green"),
                                (2/cmax, "hzdr-yellow"),
                                (3/cmax, "hzdr-orange"),
                                (4/cmax, "hzdr-red"),
                                (5/cmax, "hzdr-purple")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-rainbow2", 5.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-green"),
                                (1/cmax, "hzdr-yellow"),
                                (2/cmax, "hzdr-orange"),
                                (3/cmax, "hzdr-red"),
                                (4/cmax, "hzdr-purple"),
                                (5/cmax, "hzdr-blue")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-rainbow3", 5.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-yellow"),
                                (1/cmax, "hzdr-orange"),
                                (2/cmax, "hzdr-red"),
                                (3/cmax, "hzdr-purple"),
                                (4/cmax, "hzdr-blue"),
                                (5/cmax, "hzdr-green")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-rainbow4", 5.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-orange"),
                                (1/cmax, "hzdr-red"),
                                (2/cmax, "hzdr-purple"),
                                (3/cmax, "hzdr-blue"),
                                (4/cmax, "hzdr-green"),
                                (5/cmax, "hzdr-yellow")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn, cmax = "hzdr-rainbow5", 5.0
    d[cn] = lscm.from_list(cn, [(0/cmax, "hzdr-red"),
                                (1/cmax, "hzdr-purple"),
                                (2/cmax, "hzdr-blue"),
                                (3/cmax, "hzdr-green"),
                                (4/cmax, "hzdr-yellow"),
                                (5/cmax, "hzdr-orange")])
    d[cn+"_r"] = lscm(cn+"_r", rcms(d[cn]._segmentdata))

    cn = "hzdr-rainbow"
    d[cn] = d["hzdr-rainbow2"]
    d[cn+"_r"] =  lscm(cn+"_r", rcms(d[cn]._segmentdata))

    # Update default color map
    mpl.rcParams["image.cmap"] = "hzdr-rainbow"

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

