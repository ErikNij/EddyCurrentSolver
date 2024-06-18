#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Extended inverse distance interpolation script
# August 2013
# Dr. Thomas Wondrak (t.wondrak@hzdr.de)[basic program idea/structure]
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

from numpy import *

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def calcField(pos):
    x = pos[0] - 0.5
    y = pos[1] - 0.5
    z = pos[2] + 0.2

    mx = 0
    my = 0
    mz = 3

    v = zeros(3)

    r2 = x**2 + y**2 + z**2
    r5 = (x**2 + y**2 + z**2)**(5.0/2.0)

    s = x * mx + y * my + z * mz

    if (not (r5 == 0)):
        v[0] = (3 * x*s - mx * r2) / r5
        v[1] = (3 * y*s - my * r2) / r5
        v[2] = (3 * z*s - mz * r2) / r5

    return v




# Orginal Daten
N = 10
M = 20



# Orginal Daten
original   = zeros((N**3, 6))
cnt = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            original[cnt, 0] = 1.0 / float(N) * i
            original[cnt, 1] = 1.0 / float(N) * j
            original[cnt, 2] = 1.0 / float(N) * k

            b = calcField(original[cnt, 0:3])
            original[cnt, 3:6] = b
            cnt += 1


# Neue Punkte zur interpolation
newPoints = zeros((M**3, 3))
cnt = 0
for i in range(M):
    for j in range(M):
        for k in range(M):
            newPoints[cnt, 0] = 1.0 / float(M) * i
            newPoints[cnt, 1] = 1.0 / float(M) * j
            newPoints[cnt, 2] = 1.0 / float(M) * k
            cnt += 1


savetxt("source.dat",original)
savetxt("tnodes.dat",newPoints)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

