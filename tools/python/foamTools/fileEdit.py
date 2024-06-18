#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# File modification tools
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

import fileinput
import re

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def removeHeaderByPattern(fn, pe):

  # fn : file name
  # pe : pattern which indicates the end of the header (including pe)

  # Reset
  headEnd = False

  # Process current file
  for l in fileinput.input(fn, inplace=True):

    # Do not write lines until end of header was found
    if headEnd: sys.stdout.write (l)

    # Search end of header and remember if found
    if re.search(pe,l): headEnd = True



def removeHeaderByPrefix(fn, pf):

  # fn : file name
  # pf : prefix which indicate header lines

  # Reset
  headEnd = False

  # Process current file
  for l in fileinput.input(fn, inplace=True):

    # Write only lines without prefix
    if not re.match(pf,l): sys.stdout.write (l)



def removeBodyByPatterns(fn, ps, pe):

  # fn  : file name
  # ps : pattern which indicates the start of the body (including ps)
  # pe : pattern which indicates the end of the body (excluding pe)

  # Reset
  bodyStart = False
  bodyEnd = False

  # Process current file
  for l in fileinput.input(fn, inplace=True):

    # Head
    if not bodyStart:

      # Search start of body and do not keep line with ps
      if re.search(ps,l):
        bodyStart = True
      else:
        # Write lines
        sys.stdout.write (l)

    # Body
    elif bodyStart and not bodyEnd:

      # Search end of body and keep line with pe
      if re.search(pe,l):
        bodyEnd = True
        sys.stdout.write (l)

    # Tail
    else:

      # Write lines
      sys.stdout.write (l)



def replacePatternBySearch(fn, po, pn):

  # fn : file name
  # po : old pattern
  # pn : new pattern

  # Process current file
  for lo in fileinput.input(fn, inplace=True):

    # Search and replace patch type
    if re.search(po,lo):
      ln = lo.replace(po,pn)
    else:
      ln = lo

    # Write line
    sys.stdout.write(ln)



def replacePatternByMatch(fn, po, pn):

  # fn : file name
  # po : old pattern
  # pn : new pattern

  # Process current file
  for lo in fileinput.input(fn, inplace=True):

    # Search and replace patch type
    if re.match(po,lo):
      ln = lo.replace(po,pn)
    else:
      ln = lo

    # Write line
    sys.stdout.write(ln)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

