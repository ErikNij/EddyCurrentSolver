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

import sys
import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

if len(sys.argv) > 1:
    t=float(sys.argv[1])
else:
    print("usage: python GaInSn.py <Temperature/Celcius>")
    sys.exit(1)

# melt point
tm=10.55

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

rho=1000.0*(6.36-7.88e-4*(t-tm))
beta=0.788/rho
eta=(2.45715-0.015033*t+6.60807e-5*t*t)*0.001
nu=eta/rho
l=23.4+0.0614*(t-tm)-4.9e-5*(t-tm)**2.0
cp=368.01-0.11*t+6.67e-6*t*t
pr=eta*cp/l
s=3.33e6-3890.0*(t-tm)-48.5*(t-tm)**2.0

print("T       = %g C" % (t))
print("rho     = %g Kg/m^3" % (rho))
print("beta    = %g 1/K" % (beta))
print("eta     = %g Pa s" % (eta))
print("nu      = %g m^2/s" % (nu))
print("sigma   = %g S/m" % (s))
print("lambda  = %g W/(m K) " % (l))
print("alpha   = %g m^2/s" % (l/rho/cp))
print("Pr      = %g" % (pr))
print("cp      = %g J/kg/K" % (cp))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
