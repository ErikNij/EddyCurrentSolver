#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# March 2015
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

from numpy import array
from scipy import loadtxt

from mayavi.scripts import mayavi2
from tvtk.api import tvtk
from mayavi.modules.outline import Outline
from mayavi.modules.glyph import Glyph
from mayavi.filters.mask_points import MaskPoints
from mayavi.sources.vtk_data_source import VTKDataSource

import argparse

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# --- Argument parsing ------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Description
parser = argparse.ArgumentParser(description="This python script simply plots a vector data file with mayavi! The file must have the following structure: x,y,z,f_x,f_y,f_z,|f|")

# Optional arguments
parser.add_argument("-c", "--scale", metavar="SCALE", type=float, default=0.01, help="vector length scale factor (default: 0.05)")
parser.add_argument("-f", "--filter", metavar="FILTER", type=int, default=1, help="filter every s points while plotting (default: 2)")
parser.add_argument("-s", "--skip", metavar="SKIP", type=int, default=0, help="skip line count for plot file (default: 9)")

# Positional arguments
parser.add_argument("dataFile", metavar="DATA-FILE", nargs="?", default="plot.dat", help="Source data file for plot (default: \"plot.dat\").")

# Commence parsing
cargs = parser.parse_args()

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Load file
M = loadtxt(cargs.dataFile,skiprows=cargs.skip)

# Define column meanings
points  = M[:,0:3]
vectors = M[:,3:6]
scalars = M[:,6]

# Create and process new VTK data source
mesh = tvtk.PolyData(points=points)
mesh.point_data.vectors = vectors
mesh.point_data.scalars = scalars

@mayavi2.standalone
def te():

  # Initialize new mayavi szene
  mayavi.new_scene()

  # Load and add VTK data
  src = VTKDataSource(data = mesh)
  mayavi.add_source(src)

  # Mask filter
  m = MaskPoints()
  mayavi.add_filter(m)

  # Mask filter settings
  m.filter.set(on_ratio=cargs.filter, random_mode=True)

  # Outline module
  o = Outline()
  mayavi.add_module(o)

  # Glyph module
  g = Glyph()
  mayavi.add_module(g)

  # Glyph module settings
  g.glyph.scale_mode = "scale_by_vector"

  g.glyph.glyph.scale_factor = 0.04
  g.glyph.glyph.color_mode = "color_by_scalar"

  gs = g.glyph.glyph_source
  gs.glyph_position = "tail"
  gs.glyph_source = gs.glyph_dict["arrow_source"]

  # Readjust mask filter for mask points
  m.filter.set(random_mode=False)
  m.filter.set(random_mode=True)

te()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
