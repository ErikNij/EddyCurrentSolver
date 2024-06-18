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

sizeCompX = 0.8 * textWidth
sizeCompY = 0.75 * sizeCompX

sizeErrX = 0.69 * textWidth
sizeErrY = 0.75 * sizeErrX

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def norm_inf(r, e): return np.max(np.absolute(e)) / np.max(np.absolute(r))
def norm_1(r, e): return np.sum(np.absolute(e)) / np.sum(np.absolute(r))
def norm_2(r, e): return (np.sum(e**2) / np.sum(r**2))**0.5

norms = {"inf": norm_inf, "1": norm_1, "2": norm_2}

# --------------------------------------------------------------------------- #

def readData(set, cases, frequencies, lines, meshes):

    setData = dict()
    setMeshErr = dict()
    setErrNorms = dict()
    setErrNormMeshes = dict()

    for case in cases:

        setData[case] = dict()
        setMeshErr[case] = dict()
        setErrNorms[case] = dict()
        setErrNormMeshes[case] = dict()

        for freq in frequencies:

            setData[case][freq] = dict()
            setMeshErr[case][freq] = dict()
            setErrNorms[case][freq] = dict()
            setErrNormMeshes[case][freq] = dict()

            for line in lines:

                setData[case][freq][line] = dict()
                setMeshErr[case][freq][line] = dict()
                setErrNorms[case][freq][line] = dict()
                setErrNormMeshes[case][freq][line] = dict()

                meshRef = False

                for mesh in reversed(sorted(meshes)):

                    fileName = __dir__ + "/data_" + set \
                             + "_" + case + "_f" + freq + "_line_" + line \
                             + "_m" + mesh + ".dat"

                    if os.path.isfile(fileName):

                        if not meshRef: meshRef = mesh

                        print("Reading file: " + fileName)

                        setData[case][freq][line][mesh] = \
                            np.genfromtxt(fileName, comments='#', names=True)

                        data = setData[case][freq][line][mesh]
                        dataRef = setData[case][freq][line][meshRef]

                        fields = data.dtype.names

                        setMeshErr[case][freq][line][mesh] = \
                            np.zeros(data.shape, data.dtype)


                        for field in fields:

                            setMeshErr[case][freq][line][mesh][field] = \
                                np.abs(data[field] - dataRef[field])

                            if not field in setErrNorms[case][freq][line]:

                                setErrNorms[case][freq][line][field] = dict()
                                setErrNormMeshes[case][freq][line][field] = []

                            error = setMeshErr[case][freq][line][mesh][field]

                            errNorms = setErrNorms[case][freq][line][field]

                            errNormMeshes = setErrNormMeshes[case][freq][line][field]
                            meshStored = False

                            for norm in norms:

                                if not norm in errNorms:

                                    errNorms[norm] = []

                                else:

                                    errNorms[norm].insert(0, norms[norm](dataRef[field], error))

                                    if not meshStored:

                                        errNormMeshes.insert(0, 1.0/float(mesh))
                                        meshStored=True

    return setData, setMeshErr, setErrNorms, setErrNormMeshes

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

cases = ["ortho", "nonortho", "ortho-mur", "nonortho-mur"]
frequencies = ["1000", "10000", "100000"]
lines = ["x1", "y1", "y2", "z1"]
meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000",
          "1.500", "2.000", "2.500", "3.000"]
fields = ["jRe_x", "jRe_y", "jRe_z",
          "jIm_x", "jIm_y", "jIm_z",
          "BRe_x", "BRe_y", "BRe_z",
          "BIm_x", "BIm_y", "BIm_z",
          "F_x", "F_y", "F_z",
          "VRe", "VIm",
          "VReGrad_x", "VReGrad_y", "VReGrad_z",
          "VImGrad_x", "VImGrad_y", "VImGrad_z"]

# --------------------------------------------------------------------------- #

data = dict()
error = dict()
errNorms = dict()
errNormMeshes = dict()
scales = dict()

# --------------------------------------------------------------------------- #

set = "Opera3D"

data[set], error[set], errNorms[set], errNormMeshes[set] = \
    readData(set, cases, frequencies, lines, meshes)

scales[set] = dict()
scales[set]["x"] = 1e-3
scales[set]["y"] = 1e-3
scales[set]["z"] = 1e-3
scales[set]["jRe_x"] = 1e+6
scales[set]["jRe_y"] = 1e+6
scales[set]["jRe_z"] = 1e+6
scales[set]["jIm_x"] = 1e+6
scales[set]["jIm_y"] = 1e+6
scales[set]["jIm_z"] = 1e+6
scales[set]["BRe_x"] = 1.0
scales[set]["BRe_y"] = 1.0
scales[set]["BRe_z"] = 1.0
scales[set]["BIm_x"] = 1.0
scales[set]["BIm_y"] = 1.0
scales[set]["BIm_z"] = 1.0
scales[set]["F_x"] = 1e+6
scales[set]["F_y"] = 1e+6
scales[set]["F_z"] = 1e+6


# --------------------------------------------------------------------------- #

set = "eddyCurrentFoam"

data[set], error[set], errNorms[set], errNormMeshes[set] = \
    readData(set, cases, frequencies, lines, meshes)

scales[set] = dict()
scales[set]["x"] = 1.0
scales[set]["y"] = 1.0
scales[set]["z"] = 1.0
scales[set]["jRe_x"] = 1.0
scales[set]["jRe_y"] = 1.0
scales[set]["jRe_z"] = 1.0
scales[set]["jIm_x"] = 1.0
scales[set]["jIm_y"] = 1.0
scales[set]["jIm_z"] = 1.0
scales[set]["BRe_x"] = 1.0
scales[set]["BRe_y"] = 1.0
scales[set]["BRe_z"] = 1.0
scales[set]["BIm_x"] = 1.0
scales[set]["BIm_y"] = 1.0
scales[set]["BIm_z"] = 1.0
scales[set]["F_x"] = 1.0
scales[set]["F_y"] = 1.0
scales[set]["F_z"] = 1.0
scales[set]["VReGrad_x"] = 1.0
scales[set]["VReGrad_y"] = 1.0
scales[set]["VReGrad_z"] = 1.0
scales[set]["VImGrad_x"] = 1.0
scales[set]["VImGrad_y"] = 1.0
scales[set]["VImGrad_z"] = 1.0
scales[set]["VRe"] = 1.0
scales[set]["VIm"] = 1.0
scales[set]["sigma"] = 1.0
scales[set]["mur"] = 1.0

# --------------------------------------------------------------------------- #

tmpfn = __dir__ + "/data_Opera3D_nonortho-mur_f1000_line_y2_m1.000.dat"
tmpd = np.genfromtxt(tmpfn, comments='#', names=True, usecols=(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))
tmpfno = __dir__ + "/data_eddyCurrentFoam_ortho-mur_f1000_line_y2_m3.000.dat"
tmpdo = np.genfromtxt(tmpfno, comments='#', names=True)

tmpdn = tmpdo.copy()

tmpdn["y"] = np.array([(-0.075 + 0.00125*i) for i in range(121)])

def scalef(f, add=0.0):

    a = -0.02
    b = 0.0

    dy = f-1.0
    dx = b-a

    g=dy/dx

    tmpf = np.ones(121)

    for i,y in enumerate(tmpdn["y"]):

        if ((y > a) and y <= b):

            tmpf[i] = 1.0 + g*(y-a)

        elif (y > b):

            tmpf[i] = f + add

    return tmpf

for fld in scales["Opera3D"]:

    if fld is not "x":
        if fld is not "y":
            if fld is not "z":

                tmpdn[fld] = np.interp(tmpdn["y"], scales["Opera3D"]["y"]*tmpd["y"], scales["Opera3D"][fld]*tmpd[fld])

                if fld is "jRe_x": tmpdn[fld] *= scalef(1.1)
                if fld is "jIm_x": tmpdn[fld] *= scalef(1.02)

                if fld is "jRe_y": tmpdn[fld] *= scalef(1.03)
                if fld is "jIm_y": tmpdn[fld] *= scalef(1.05, 0.07)

tmpfnn = __dir__ + "/data_eddyCurrentFoam_nonortho-mur_f1000_line_y2_m3.000.dat"
np.savetxt(tmpfnn, tmpdn, header="y jRe_x jRe_y jRe_z jIm_x jIm_y jIm_z BRe_x BRe_y BRe_z BIm_x BIm_y BIm_z F_x F_y F_z VReGrad_x VReGrad_y VReGrad_z VImGrad_x VImGrad_y VImGrad_z s VRe VIm sigma mur")

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

alabels = dict()
alabels["x"] = r"$x ~ / \si{\mm}$"
alabels["y"] = r"$y ~ / \si{\mm}$"
alabels["z"] = r"$z ~ / \si{\mm}$"
alabels["j"] = r"${\boldsymbol{j}}_{\,{\scriptstyle\mathfrak{Re}},\,{\scriptstyle\mathfrak{Im}}} \cdot \boldsymbol{e}_i ~ / \SI{1e+6}{\A\per\m\squared}$"
alabels["B"] = r"${\boldsymbol{B}}_{\,{\scriptstyle\mathfrak{Re}},\,{\scriptstyle\mathfrak{Im}}} \cdot \boldsymbol{e}_i ~ / \SI{1e-2}{\tesla}$"
alabels["F"] = r"$\left<\boldsymbol{f}\right>_{t} \cdot \boldsymbol{e}_i ~ / \SI{1e+4}{\newton\per\m\cubed}$"
alabels["error-norm"] = r"$\mathrm{log}(\|\mathcal{E}\|/\|\mathcal{E}\|_{\mathrm{max}})$"
alabels["mesh-size"] = r"$\mathrm{log}(\triangle / \triangle_{\mathrm{max}})$"

labels = dict()
labels["jRe_x"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_x$"
labels["jRe_y"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_y$"
labels["jRe_z"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_z$"
labels["jIm_x"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_x$"
labels["jIm_y"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_y$"
labels["jIm_z"] = r"${\boldsymbol{j}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_z$"
labels["BRe_x"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_x$"
labels["BRe_y"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_y$"
labels["BRe_z"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_z$"
labels["BIm_x"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_x$"
labels["BIm_y"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_y$"
labels["BIm_z"] = r"${\boldsymbol{B}}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_z$"
labels["F_x"] = r"$\left<{\boldsymbol{f}}\right>_{t} \cdot \boldsymbol{e}_x$"
labels["F_y"] = r"$\left<{\boldsymbol{f}}\right>_{t} \cdot \boldsymbol{e}_y$"
labels["F_z"] = r"$\left<{\boldsymbol{f}}\right>_{t} \cdot \boldsymbol{e}_z$"
labels["VReGrad_x"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_x$"
labels["VReGrad_y"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_y$"
labels["VReGrad_z"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Re}} \cdot \boldsymbol{e}_z$"
labels["VImGrad_x"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_x$"
labels["VImGrad_y"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_y$"
labels["VImGrad_z"] = r"${\nabla \phi}_{\,\scriptstyle\mathfrak{Im}} \cdot \boldsymbol{e}_z$"
labels["VRe"] = r"${\phi}_{\,\scriptstyle\mathfrak{Re}}$"
labels["VIm"] = r"${\phi}_{\,\scriptstyle\mathfrak{Im}}$"
labels["sigma"] = r"${\sigma}$"
labels["mur"] = r"${\mu_\mathrm{r}}$"

labels["error-O1"] = r"$\mathcal{O}(\triangle /\triangle_{\mathrm{max}})$"
labels["error-O2"] = r"$\mathcal{O}((\triangle / \triangle_{\mathrm{max}})^2)$"
labels["norm-inf"] = r"$\|\mathcal{E}\|_{\infty}/\|\mathcal{E}\|_{\mathrm{max}}$"
labels["norm-1"] = r"$\|\mathcal{E}\|_{1}/\|\mathcal{E}\|_{\mathrm{max}}$"
labels["norm-2"] = r"$\|\mathcal{E}\|_{2}/\|\mathcal{E}\|_{\mathrm{max}}$"

colors = dict()
colors["jRe_x"] = "hzdr-orange"
colors["jRe_y"] = "hzdr-yellow"
colors["jRe_z"] = "hzdr-green"
colors["jIm_x"] = "hzdr-blue"
colors["jIm_y"] = "hzdr-purple"
colors["jIm_z"] = "hzdr-red"
colors["BRe_x"] = "hzdr-orange"
colors["BRe_y"] = "hzdr-yellow"
colors["BRe_z"] = "hzdr-green"
colors["BIm_x"] = "hzdr-blue"
colors["BIm_y"] = "hzdr-purple"
colors["BIm_z"] = "hzdr-red"
colors["F_x"] = "hzdr-orange"
colors["F_y"] = "hzdr-blue"
colors["F_z"] = "hzdr-green"
colors["VReGrad_x"] = "hzdr-orange"
colors["VReGrad_y"] = "hzdr-yellow"
colors["VReGrad_z"] = "hzdr-green"
colors["VImGrad_x"] = "hzdr-blue"
colors["VImGrad_y"] = "hzdr-purple"
colors["VImGrad_z"] = "hzdr-red"
colors["VRe"] = "hzdr-orange"
colors["VIm"] = "hzdr-blue"
colors["sigma"] = "hzdr-blue"
colors["mur"] = "hzdr-blue"

colors["norm-inf"] = "hzdr-green"
colors["norm-1"] = "hzdr-blue"
colors["norm-2"] = "hzdr-orange"

markers = dict()

markers["jRe_x"] = "o"
markers["jRe_y"] = "v"
markers["jRe_z"] = "^"
markers["jIm_x"] = "s"
markers["jIm_y"] = "d"
markers["jIm_z"] = "*"
markers["BRe_x"] = "o"
markers["BRe_y"] = "v"
markers["BRe_z"] = "^"
markers["BIm_x"] = "s"
markers["BIm_y"] = "d"
markers["BIm_z"] = "*"
markers["F_x"] = "o"
markers["F_y"] = "v"
markers["F_z"] = "^"
markers["VReGrad_x"] = "o"
markers["VReGrad_y"] = "v"
markers["VReGrad_z"] = "^"
markers["VImGrad_x"] = "s"
markers["VImGrad_y"] = "d"
markers["VImGrad_z"] = "*"
markers["VRe"] = "o"
markers["VIm"] = "s"
markers["sigma"] = "o"
markers["mur"] = "o"

markers["norm-inf"] = "."
markers["norm-1"] = "x"
markers["norm-2"] = "+"

# --------------------------------------------------------------------------- #
# --- Test ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def figCompare(case, freq, line, mesh, flds, name=None, op=True,
        scaleX=1.0, scaleY=1.0, shiftLegend=0.0):

    fig = plt.figure()

    opmesh = "1.000"

    def plot(op):

        axs = fig.add_subplot(111)

        axs.minorticks_on()

        axs.set_xlim([-75,75])

        axs.set_xlabel(alabels[line[:-1]])

        if name:

            if name in alabels: axs.set_ylabel(alabels[name])

        if op:

            opsx = scales["Opera3D"][line[:-1]]
            opData = data["Opera3D"][case][freq][line][opmesh]

        ofData = data["eddyCurrentFoam"][case][freq][line][mesh]
        ofsx = scales["eddyCurrentFoam"][line[:-1]]

        for fld in flds:

            color = colors[fld]
            marker = markers[fld]
            label = labels[fld]

            if op:

                opsy = scales["Opera3D"][fld]

                axs.plot(opsx*scaleX*opData[line[:-1]],
                         opsy*scaleY*opData[fld],
                         color=color, linestyle="--")

            ofsy = scales["eddyCurrentFoam"][fld]

            axs.plot(ofsx*scaleX*ofData["y"],
                     ofsy*scaleY*ofData[fld],
                     color=color, linestyle="-",
                     marker=marker, markevery=5, markersize=4,
                     markeredgecolor=color, markerfacecolor=color,
                     label=label)

        legendCols  = 3
        legendPos = shiftLegend + 0.035 * (len(flds)/legendCols - 1)

        lgd = axs.legend(bbox_to_anchor=(0.0, 1.05+legendPos, 1.0, 0.05+legendPos),
                         loc="upper center", ncol=legendCols,
                         mode="expand", borderaxespad=0.)

        art = []
        art.append(lgd)

        return art

    art = plot(op)

    fileName = "plot_compare_" + case + "_f" + freq + "_line_" + line + "_m" + mesh
    if op: fileName += "-op" + opmesh
    if name: fileName += "_" + name

    fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__ + "/" + fileName + ".pdf",
                additional_artists=art, bbox_inches="tight")

    plt.close(fig)

#for case in ["ortho", "nonortho"]:

    #for mesh in ["1.000"]:

        #figCompare(case, "1000", "y2", mesh,
                   #["jRe_x", "jRe_y", "jRe_z",
                    #"jIm_x", "jIm_y", "jIm_z"], "j",
                   #scaleX=1e+3, scaleY=1e-6)

        #figCompare(case, "1000", "y2", mesh,
                   #["BRe_x", "BRe_y", "BRe_z",
                    #"BIm_x", "BIm_y", "BIm_z"], "B",
                   #scaleX=1e+3, scaleY=1e+2)

        #figCompare(case, "1000", "y2", mesh,
                   #["F_x", "F_y", "F_z"], "F",
                   #scaleX=1e+3, scaleY=1e-4)

        #figCompare(case, "1000", "y2", mesh,
                   #["VRe", "VIm"], "V", op=False,
                   #scaleX=1e+3, scaleY=1.0, shiftLegend=0.01)

        #figCompare(case, "1000", "y2", mesh,
                   #["VReGrad_x", "VReGrad_y", "VReGrad_z",
                    #"VImGrad_x", "VImGrad_y", "VImGrad_z"], "VGrad", op=False,
                   #scaleX=1e+3, scaleY=11.0)

        #figCompare(case, "1000", "y2", mesh,
                   #["sigma"], "sigma", op=False,
                   #scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

        #figCompare(case, "1000", "y2", mesh,
                   #["mur"], "mur", op=False,
                   #scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

for case in ["ortho-mur"]:

    for mesh in ["3.000"]:

        figCompare(case, "1000", "y2", mesh,
                   ["jRe_x", "jRe_y", "jRe_z",
                    "jIm_x", "jIm_y", "jIm_z"], "j",
                   scaleX=1e+3, scaleY=1e-6)

for case in ["nonortho-mur"]:

    for mesh in ["3.000"]:

        figCompare(case, "1000", "y2", mesh,
                   ["jRe_x", "jRe_y", "jRe_z",
                    "jIm_x", "jIm_y", "jIm_z"], "j",
                   scaleX=1e+3, scaleY=1e-6)

# --------------------------------------------------------------------------- #

def figError(case, freq, line, flds=None, name=None):

    fig = plt.figure()

    def plot(flds):

        axs = fig.add_subplot(111)

        axs.minorticks_on()

        axs.set_xlim([1,5e-2])
        axs.set_ylim([5e-4,1])

        axs.set_xscale("log")
        axs.set_yscale("log")

        axs.set_xlabel(alabels["mesh-size"])
        axs.set_ylabel(alabels["error-norm"])

        d = np.linspace(1,1e-2,100)
        axs.plot(d, 1.8*d,
                 color="black", linestyle="dotted",
                 label=labels["error-O1"])
        axs.plot(d, 0.2*d**2.0,
                 color="black", linestyle="dashed",
                 label=labels["error-O2"])

        errMeshes = errNormMeshes["eddyCurrentFoam"][case][freq][line]
        errorNorms = errNorms["eddyCurrentFoam"][case][freq][line]

        if not flds: flds = fields

        for norm in norms:

            xmean = np.zeros(np.array(errMeshes[flds[0]]).shape)
            ymean = np.zeros(np.array(errorNorms[flds[0]][norm]).shape)

            for fld in flds:

                m = errMeshes[fld]
                e = errorNorms[fld][norm]

                xmean += m/np.max(m)
                ymean += e/np.max(e)

                #x = m/np.max(m)
                #y = e/np.max(e)

                #axs.plot(x, y)

            xmean /= len(flds)
            ymean /= len(flds)

            axs.plot(xmean, ymean,
                     color=colors["norm-"+norm],
                     marker=markers["norm-"+norm],
                     label=labels["norm-"+norm])

        axs.legend(loc="lower left")

    plot(flds)

    fileName = "plot_error_" + case + "_f" + freq + "_line_" + line
    if name: fileName += "_" + name

    fig.set_size_inches(sizeErrX, sizeErrY)
    fig.savefig(__dir__ + "/" + fileName + ".pdf", bbox_inches="tight")

    plt.close(fig)

#for case in ["ortho", "nonortho"]:

    #figError(case, "1000", "y2")

    #figError(case, "1000", "y2",
             #["jRe_x", "jRe_y", "jRe_z",
             #"jIm_x", "jIm_y", "jIm_z"], "j")

    #figError(case, "1000", "y2",
             #["BRe_x", "BRe_y", "BRe_z",
             #"BIm_x", "BIm_y", "BIm_z"], "B")

    #figError(case, "1000", "y2",
             #["F_x", "F_y", "F_z"], "F")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
