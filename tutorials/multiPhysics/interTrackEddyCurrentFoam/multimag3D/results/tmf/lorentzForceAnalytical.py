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

import math as m
import numpy as np

from foamTools.ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

k=2.0*m.pi/0.4476
H0 = 0.03
R0 = 0.03
nr = 61
nz = 121

mu=4.0e-7*m.pi

rho=6352.55
sigma = 3.28891e+06
eta=0.00218292
nu=eta/rho
omega = 2.0*m.pi*50.0
#B0 = 1.988e-03
#B0 = 1.99984e-03
B0 = 2.032e-03

print("k          : {} 1/m".format(k))
print("H0, R0     : {}, {} m".format(H0, R0))
print("nr, nz     : {},{}".format(nr, nz))

print("rho        : {} kg/m^3".format(rho))
print("sigma      : {} S/m".format(sigma))
print("eta        : {} Pa s".format(eta))
print("nu         : {} m^2/s".format(nu))
print("omega      : {} 1/s".format(omega))
print("B0         : {} T".format(B0))

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Shielding parameter
def S():

    return mu*sigma*omega*m.pow(R0,2)



# Forcing parameter
def F():

    return sigma*omega*k*m.pow(B0,2)*m.pow(R0,5)/(4.0*rho*m.pow(nu,2))



# Lorentz-Force density
def fl(r,z):

    return 0.125*sigma*omega*k*m.pow(B0,2)*m.pow(r,2)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Dimless parameters
print("S          : {} ".format(S()))
print("F          : {} ".format(F()))

# Mesh and discretisation
rl = np.linspace(0.0,R0,nr)
zl = np.linspace(-H0,H0,nz)

R,Z = np.meshgrid(rl,zl,indexing="ij")
F   = np.zeros(R.shape)


# Calculate and write data
with open(__dir__+"/lorentzForceAnalytical.dat","w") as lff:

    lff.write("# Variables: " + "r" + " " + "z" + " " + "F" + "\n")

    lff.write("# Resolution: " + "nr=" + str(nr) + ", " + "nz=" + str(nz) + "\n")

    for ri in range(nr):
        for zi in range(nz):

            F[ri,zi] = fl(rl[ri],zl[zi])

            lff.write(str(rl[ri]) + " " + str(zl[zi]) + " " + str(F[ri,zi]) + "\n")

print("Fmin, Fmax : {},{}".format(np.min(F), np.max(F)))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
