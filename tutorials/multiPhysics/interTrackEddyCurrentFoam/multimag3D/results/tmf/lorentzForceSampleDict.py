#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
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

import numpy as np

from foamTools.ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

H0 = 0.029999
R0 = 0.029999

nr = 61
nz = 121

print("H0, R0     : {}, {} m".format(H0, R0))
print("nr, nz     : {},{}".format(nr, nz))

# Mesh and discretisation
rl = np.linspace(0.0,R0,nr)
zl = np.linspace(-H0,H0,nz)

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    with open(__dir__+"/lorentzForce.samleDict", "w") as f:

        # Define short indented line with line break
        def i(iL,cS,eS="\n"): return objectIndent(cS + eS,iLevel=iL)

        f.write(objectHeader("sampleDict", "dictionary"))

        f.write(i(0, "interpolationScheme cellPointFace;"))
        f.write(i(0, "setFormat           raw;"))
        f.write(i(0, "surfaceFormat       null;"))

        f.write(i(0, "sets "))
        f.write(i(0, "("))
        f.write(i(1, "planexz"))
        f.write(i(1, "{"))
        f.write(i(2, "type cloud;"))
        f.write(i(2, "axis    xyz;"))
        f.write(i(2, "points"))
        f.write(i(2, "("))

        z0 = 0.03

        for ri in range(len(rl)):
            for zi in range(len(zl)):

                f.write(i(3,"(" + str(rl[ri]) + " 0.0 " + str(z0 + zl[zi]) + ")"))

        f.write(i(2, ");"))
        f.write(i(1, "}"))
        f.write(i(0, ");"))

        f.write("\n")

        f.write(i(0, "surfaces "))
        f.write(i(0, "("))
        f.write(i(0, ");"))

        f.write("\n")

        f.write(i(0, "fields "))
        f.write(i(0, "("))
        f.write(i(1, "F"))
        f.write(i(0, ");"))

        f.write("\n")

        f.write(objectFooter())

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
