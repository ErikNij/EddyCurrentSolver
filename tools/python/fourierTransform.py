#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Compute the discrete fourier transform of data from file
# Jan 2018
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

import numpy as np
import pylab as p

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

print(sys.argv)

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    print("Missing file name.")
    sys.exit(1)

if len(sys.argv) > 2:
    deltaT = float(sys.argv[2])
else:
    print("Missing time step width.")
    sys.exit(1)

if len(sys.argv) > 3:
    omitlT = float(sys.argv[3])
else:
    print("Missing omitl time.")
    sys.exit(1)

if len(sys.argv) > 4:
    omitrT = float(sys.argv[4])
else:
    print("Missing omitr time.")
    sys.exit(1)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

with open(filename) as file:
     floats = map(float,file)

# data
UX = np.array(floats)
# number of data values
NT = len(floats)
NTb2 = int(NT/2)

# mean
UXmean = UX.mean()
print("data average       : %g" % (UXmean))
UX = UX - UXmean

# time step DT = (tmax-tmin)/(NT-1)
DT = deltaT
print("delta_t            : %g" % (DT))

# smallest frequency
DF = 1.0/((NT-1.0)*DT)
print("smallest frequency : %g" % (DF))

# real fourier transform
fft=(2.0/NT)*np.fft.rfft(UX)

# power spectrum
power=abs(fft)

# discrete freqs
if NT%2==0:
   F = np.linspace(0,NT*DF/2,NTb2+1)
else:
   F = np.linspace(0,(NT-1)*DF/2,(NT-1)/2+1)

# omit left crap
omitlFreq = 1.0/omitlT
omitlPeriode = np.float64(1.0) / omitlFreq
omitl = 0
while F[omitl] < omitlFreq and omitl < NTb2:
    omitl += 1
print("omitl              : %d" % (omitl))
print("omitlFreq          : %g Hz" % (omitlFreq))
print("omitlPeriode       : %g s" % (omitlPeriode))

# omit right crap
omitrFreq = 1.0/omitrT
omitrPeriode = np.float64(1.0) / omitrFreq
omitr = 0
while F[omitr] < omitrFreq and omitr < NTb2:
    omitr += 1
print("omitr              : %d" % (omitr))
print("omitrFreq          : %g Hz" % (omitrFreq))
print("omitrPeriode       : %g s" % (omitrPeriode))

# find max
maxi = np.argmax(power[omitl:omitr]) + omitl
maxiFreq = F[maxi]
maxiPeriode = np.float64(1.0) / maxiFreq
print("maxi               : %d" % (maxi))
print("maxiFreq           : %g Hz" % (maxiFreq))
print("maxiPeriode        : %g s" % (maxiPeriode))

# cut small freqs
cutFactor = 0.1
cutWindow = 10
cut = NTb2
while power[cut-cutWindow:cut].mean() < cutFactor*power[maxi] \
    and cut-cutWindow > maxi:
    cut -= 1
cutFreq = F[cut]
cutPeriode = np.float64(1.0) / cutFreq
print("cut                : %d" % (cut))
print("cutFreq            : %g Hz" % (cutFreq))
print("cutPeriode         : %g s" % (cutPeriode))
print("cutFactor          : %g" % (cutFactor))
print("cutWindow          : %d" % (cutWindow))

print("N                  : %d" % (NT))
print("Peak frequency, f  : %g Hz" % (maxiFreq))
print("Periode, T         : %g s"  % (1.0/maxiFreq))

p.figure(figsize=(8,10))

ax1 = p.subplot(211)
p.title(filename)
p.xlabel("t / s")
p.ylabel("signal")
p.grid()
p.plot(DT*np.linspace(0,NT,NT), UX+UXmean)

ax2 = p.subplot(212)
p.xlabel("f / Hz")
p.ylabel("power spectrum")
p.xlim(0,cutFreq)
p.ylim(0,1.5*power[maxi])
p.grid()
p.plot(F,power)
p.axvline(omitlFreq,c="black",ls="dashed")
p.axvline(omitrFreq,c="black",ls="dashed")
p.axvline(maxiFreq,c="red",ls="solid")
p.annotate("N = %d, DT = %g\nf = %g Hz\nT = %g s" % (NT, DT, maxiFreq, (1.0/maxiFreq)), (1.2*maxiFreq,1.15*power[maxi]), color="red")

p.savefig(os.path.dirname(filename) + "/" + os.path.splitext(os.path.basename(filename))[0] + ".pdf")

#p.show()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
