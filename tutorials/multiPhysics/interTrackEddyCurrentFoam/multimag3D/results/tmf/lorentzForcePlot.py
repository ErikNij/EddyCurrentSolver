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

import matplotlib.pyplot as plt

import plotTools.latex as latex
import plotTools.hzdr as hzdr

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

nr = 61
nz = 121

print("nr, nz     : {},{}".format(nr, nz))

# --------------------------------------------------------------------------- #

fontsize = 11
fontfamily = "serif" # only for clabels
locale="de_DE.utf8"
preamble = ["\\usepackage[utf8x]{inputenc}",
            "\\usepackage[T1]{fontenc}",
            "\\usepackage{amsmath,amssymb,amsthm,amsfonts,mathrsfs}",
            "\\usepackage{mathtools}",
            "\\usepackage{newtxtext,newtxmath}",
            "\\usepackage[locale=" + locale[3:5] + "]{siunitx}"]

latex.pdflatexify(fontsize=fontsize, fontfamily=fontfamily,
                  locale=locale, preamble=preamble)

hzdr.colors()

# --------------------------------------------------------------------------- #

textWidth = (210-25-25)/25.4

sizeCompX = 0.49 * textWidth
sizeCompY = 2.0 * sizeCompX

sizeErrX = 0.69 * textWidth
sizeErrY = 0.75 * sizeErrX

baseName = "lorentzForce"

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def norm_inf(r, e): return np.max(np.absolute(e)) / np.max(np.absolute(r))
def norm_1(r, e): return np.sum(np.absolute(e)) / np.sum(np.absolute(r))
def norm_2(r, e): return (np.sum(e**2) / np.sum(r**2))**0.5

norms = {"inf": norm_inf, "1": norm_1, "2": norm_2}

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

data = dict()

R = dict() # Radius
Z = dict() # Height

F = dict() # Force

E = dict() # Error
D = dict() # Delta of grid
N = dict() # Norm of error

# --------------------------------------------------------------------------- #

set = "Analytical"

# Read data
if True:

    fileName = __dir__+"/"+baseName+set+".dat"

    print("Reading file: " + fileName)

    data[set] = np.genfromtxt(fileName, comments="#")

    R[set]    = data[set][:,0].reshape(nr,nz)
    Z[set]    = data[set][:,1].reshape(nr,nz)

    F[set]    = [np.zeros(R[set].shape) for i in range(3)]
    F[set][2] = data[set][:,2].reshape(nr,nz)

    # Scale to mm
    R[set] = 1000.0 * R[set]
    Z[set] = 1000.0 * Z[set]

# --------------------------------------------------------------------------- #

set = "Opera3D"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["coarse", "fine"]

# Read data
for mesh in meshes:

    fileName = __dir__+"/"+baseName+set+"_"+mesh+".dat"

    print("Reading file: " + fileName)

    data[set][mesh] = np.genfromtxt(fileName, comments="#")

    R[set][mesh]    = 1000.0 * data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]    = 1000.0 * data[set][mesh][:,1].reshape(nr,nz)

    F[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]
    F[set][mesh][0] = data[set][mesh][:,2].reshape(nr,nz)
    F[set][mesh][2] = data[set][mesh][:,3].reshape(nr,nz)

# --------------------------------------------------------------------------- #

set = "EddyCurrentFoam"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000", "1.500", "2.000", "2.500"]

# Read data
for mesh in meshes:

    fileName = __dir__+"/"+baseName+set+"_"+mesh+".dat"

    print("Reading file: " + fileName)

    data[set][mesh] = np.genfromtxt(fileName, comments="#")

    R[set][mesh]    = data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]    = data[set][mesh][:,2].reshape(nr,nz)

    F[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]
    F[set][mesh][0] = data[set][mesh][:,3].reshape(nr,nz)
    F[set][mesh][1] = data[set][mesh][:,4].reshape(nr,nz)
    F[set][mesh][2] = data[set][mesh][:,5].reshape(nr,nz)

    # Scale to mm
    R[set][mesh]    = 1000.0 * (R[set][mesh])
    Z[set][mesh]    = 1000.0 * (Z[set][mesh] - 0.03)

    # Flip y-axis
    F[set][mesh][1] = -F[set][mesh][1]

# Calculate errors with repect to last mesh in list
for mesh in meshes[:-1]:

    E[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]

    for i in [0, 2]:

        E[set][mesh][i] = abs(F[set][meshes[-1]][i] - F[set][mesh][i])

# Init norms
for norm in norms.keys():

    N[set][norm] = dict()

    for i in range(3):

        N[set][norm][i] = dict()

# Calculate norms and global maximum
Dmax = 0.0
Nmax = [ 0.0 for i in range(3)]
for mesh in meshes[:-1]:

    D[set][mesh] = 1.0/float(mesh)
    Dmax = max(Dmax, D[set][mesh])

    for norm in norms.keys():

        for i in [0, 2]:

            N[set][norm][i][mesh] = norms[norm](F[set][meshes[-1]][i], E[set][mesh][i])
            Nmax[i] = max(Nmax[i], N[set][norm][i][mesh])

# Scale norms and mesh delta
for mesh in meshes[:-1]:

    D[set][mesh] /= Dmax

    for norm in norms.keys():

        for i in [0, 2]:

            N[set][norm][i][mesh] /= Nmax[i]

# --------------------------------------------------------------------------- #

set = "EddyCurrentFoam_lowf"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000", "1.500", "2.000"]

# Read data
for mesh in meshes:

    fileName = __dir__+"/"+baseName+set+"_"+mesh+".dat"

    print("Reading file: " + fileName)

    data[set][mesh] = np.genfromtxt(fileName, comments="#")

    R[set][mesh]    = data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]    = data[set][mesh][:,2].reshape(nr,nz)

    F[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]
    F[set][mesh][0] = data[set][mesh][:,3].reshape(nr,nz)
    F[set][mesh][1] = data[set][mesh][:,4].reshape(nr,nz)
    F[set][mesh][2] = data[set][mesh][:,5].reshape(nr,nz)

    # Scale to mm
    R[set][mesh]    = 1000.0 * (R[set][mesh])
    Z[set][mesh]    = 1000.0 * (Z[set][mesh] - 0.03)

    # Flip y-axis
    F[set][mesh][1] = -F[set][mesh][1]

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

labelAxisR = r"$r ~ / \si{\mm}$"
labelAxisZ = r"$z ~ / \si{\mm}$"

labelAxisE = r"$\mathrm{log}(\|\mathcal{E}\|/\|\mathcal{E}\|_{\mathrm{max}})$"
labelAxisD = r"$\mathrm{log}(\triangle / \triangle_{\mathrm{max}})$"
labelAxisS = r"$\mathrm{log}(R_\mathrm{\infty}/R_J)$"

labelE = dict()
labelE["O1"] = r"$\mathcal{O}(\triangle /\triangle_{\mathrm{max}})$"
labelE["O2"] = r"$\mathcal{O}((\triangle / \triangle_{\mathrm{max}})^2)$"
labelE["Os3"] = r"$\mathcal{O}((R_\mathrm{\infty}/R_J)^{-3})$"
labelE["OsConst"] = r"$R_\mathrm{\infty}/R_J = \mathrm{const.}$"
labelE["inf"] = r"$\|\mathcal{E}\|_{\infty}/\|\mathcal{E}\|_{\mathrm{max}}$"
labelE["1"] = r"$\|\mathcal{E}\|_{1}/\|\mathcal{E}\|_{\mathrm{max}}$"
labelE["2"] = r"$\|\mathcal{E}\|_{2}/\|\mathcal{E}\|_{\mathrm{max}}$"

markerE = dict()
markerE["O1"] = ""
markerE["O2"] = ""
markerE["Os3"] = ""
markerE["inf"] = "."
markerE["1"] = "x"
markerE["2"] = "+"

levels = [np.array([-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5]),
          np.array([0]),
          np.array([1,2,3,4,5,6,7,8])]
#levels = [np.linspace(-4.00, 0.00,9),
          #np.linspace( 0.00, 0.00,1),
          #np.linspace( 0.00, 8.00,9)]

#colors = "black"

plots = dict()

# --------------------------------------------------------------------------- #
# --- Analytical solution --------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([-30,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Analytical"
        c = ax.contour(R[set], Z[set], F[set][2],
                       levels=levels[2])
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        [l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "Fz")

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "Analytical")

# --------------------------------------------------------------------------- #
# --- Low frequency approximation ------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([-30,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Analytical"
        ax.contour(R[set], Z[set], F[set][2],
                   levels=levels[2], linestyles="dotted")

        set = "EddyCurrentFoam_lowf"
        mesh = "1.000"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][2],
                       levels=levels[2], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        [l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "Fz")

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "LowFrequencyFz")

# --------------------------------------------------------------------------- #
# --- Contour plots --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

Opera3D_ref = "fine"
EddyCurrentFoam_ref = "1.000"

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([-30,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Opera3D"
        mesh = Opera3D_ref
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][0],
                       levels=levels[0], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = EddyCurrentFoam_ref
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][0],
                       levels=levels[0], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        [l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs,"Fr")

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonFr")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([-30,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Analytical"
        ax.contour(R[set], Z[set], F[set][2],
                   levels=levels[2], linestyles="dotted")

        set = "Opera3D"
        mesh = Opera3D_ref
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][2],
                       levels=levels[2], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = EddyCurrentFoam_ref
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][2],
                       levels=levels[2], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        [l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "Fz")

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonFz")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([-30,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Opera3D"
        mesh = Opera3D_ref
        magF = (F[set][mesh][0]**2 + F[set][mesh][1]**2 + F[set][mesh][2]**2)**0.5
        c = ax.contour(R[set][mesh], Z[set][mesh], magF,
                       levels=levels[2], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = EddyCurrentFoam_ref
        magF = (F[set][mesh][0]**2 + F[set][mesh][1]**2 + F[set][mesh][2]**2)**0.5
        c = ax.contour(R[set][mesh], Z[set][mesh], magF,
                       levels=levels[2], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        [l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "F")

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonF")

# --------------------------------------------------------------------------- #
# --- Error plots ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([5e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "Fr")

    fig.set_size_inches(sizeErrX, sizeErrY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorFr")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([5e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "Fz")

    fig.set_size_inches(sizeErrX, sizeErrY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorFz")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([5e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        magnMax = 0.0
        magn = dict()
        for norm in sorted(norms.keys()):
            n0 = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            n2 = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            magn[norm] = (n0**2 + n2**2)**0.5
            for m in magn[norm]:
                magnMax = max(magnMax, m)

        for norm in sorted(norms.keys()):
            n = magn[norm]/magnMax
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "F")

    fig.set_size_inches(sizeErrX, sizeErrY)
    fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorF")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
