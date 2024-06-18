#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)
#
#  Example usage
#  ~~~~~~~~~~~~~
#
#  coil_scale      = 1e-3
#
#  coil_n          = 10
#  coil_step       = 14.9
#  coil_origin     = [0.0, 0.0, 7.0]
#
#  coil_path       = {"shape": "loop",
#                     "n":     36,
#                     "r":     105.0}
#  coil_path       = {"shape": "racetrack",
#                     "n":     9,
#                     "r":     10.0,
#                     "x":     50.0,
#                     "y":     100.0}
#
#  coil_bundle     = {"shape": "point"}
#  coil_bundle     = {"shape": "circle",
#                     "n":     36,
#                     "r":     5.0}
#  coil_bundle     = {"shape": "rectangle",
#                     "n":     10,
#                     "r":     10.0,
#                     "z":     8.0}
#
#  coil_current    = m.sqrt(2.0) * 260.0
#  coil_nNonOrto   = 10
#  coil_frequency  = 6300.0
#
#  coils = inductorCoils("ARRAY", csn, par.coil_bundle, par.coil_path,
#                        par.coil_current, par.coil_n, par.coil_step,
#                        origin=par.coil_origin, axis=2, scale=par.coil_scale)
#
#  writeCoilFeatureEdgeMeshes(par.dir_case, coils)
#  writeEdgeBiotSavartProperties(par.dir_case, coils, par.coil_nNonOrto)
#  writeFrequency(par.dir_case, par.coil_frequency)
#

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

import math as m
import numpy as np

from foamTools.ioInfo import objectIndent, objectHeader, objectFooter
from foamTools.math import rotationMatrix

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# +++ Paths +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

def edgeLoopFromPoints(points, edgeStart=0):

    edges = list()

    l = len(points) - 1

    for i in range(l):

        I = edgeStart + i

        edges.append([I, I + 1])

    edges.append([edgeStart + l, edgeStart])

    return edges



def pathLoop(pathDict, bundleDict, filamentI, edgeStart=0):
    """

    pathDict-Keys
    ----------
    n : int, Number of edges
    r : float, Coil loop radius
    """

    if not "n" in pathDict:

        raise KeyError("Number of edges (n) is missing.")

    if not isinstance(pathDict["n"], int):

        raise KeyError("Number of edges (n) needs to be of type int.")

    if not pathDict["n"] > 0:

        raise ValueError("Number of edges (n) needs to be larger than 0.")

    if not "r" in pathDict:

        raise KeyError("Coil loop radius (r) is missing.")

    if pathDict["r"] <= 0.0:

        raise ValueError("Coil loop radius (r) must be positive.")

    f = bundle[bundleDict["shape"]](bundleDict, filamentI)

    points = list()

    r = f[0] + pathDict["r"]
    z = f[1]
    phii = 1.0/pathDict["n"] * 2.0*m.pi

    for i in range(pathDict["n"]):

        p = np.zeros(3)

        p[0] = r * m.cos(i * phii)
        p[1] = r * m.sin(i * phii)
        p[2] = z

        points.append(p)

    edges = edgeLoopFromPoints(points, edgeStart)

    return points, edges



def pathRaceTrack(pathDict, bundleDict, filamentI, edgeStart=0):
    """

    pathDict-Keys
    ----------
    n : int, Number of edges for each corner arc
    r : float, Coil corner radius
    x : float, Coil size in x-direction
    y : float, Coil size in y-direction
    """

    if not "n" in pathDict:

        raise KeyError("Number of edges (n) is missing.")

    if not isinstance(pathDict["n"], int):

        raise KeyError("Number of edges (n) needs to be of type int.")

    if not pathDict["n"] > 0:

        raise ValueError("Number of edges (n) needs to be larger than 0.")

    if not "r" in pathDict:

        raise KeyError("Inner coil corner radius (r) is missing.")

    if pathDict["r"] <= 0.0:

        raise ValueError("Inner coil corner radius (r) must be positive.")

    if not "x" in pathDict:

        raise KeyError("Coil width (x) is missing.")

    if not "y" in pathDict:

        raise KeyError("Coil height (y) is missing.")

    if pathDict["x"] <= 0.0 or pathDict["y"] <= 0.0:

        raise ValueError("Coil sizes (x/y) must be positive.")

    f = bundle[bundleDict["shape"]](bundleDict, filamentI)

    points = list()

    x = pathDict["x"] - pathDict["r"]
    y = pathDict["y"] - pathDict["r"]
    r = f[0] +  pathDict["r"]
    z = f[1]
    phii = 1.0/pathDict["n"] * m.pi/2.0
    phi0 = m.pi * np.array([0.0, 0.5, 1.0, 1.5])
    s = [[1, 1], [-1, 1], [-1, -1], [1, -1]]

    for c in range(len(s)):

        for i in range(pathDict["n"] + 1):

            p = np.zeros(3)

            p[0] = s[c][0]*x + r * m.cos(phi0[c] + i*phii)
            p[1] = s[c][1]*y + r * m.sin(phi0[c] + i*phii)
            p[2] = z

            points.append(p)

    edges = edgeLoopFromPoints(points, edgeStart)

    return points, edges



def pathRectangle(pathDict, bundleDict, filamentI, edgeStart=0):
    """

    pathDict-Keys
    ----------
    x : float, Coil size in x-direction
    y : float, Coil size in y-direction
    """

    if not "x" in pathDict:

        raise KeyError("Coil width (x) is missing.")

    if not "y" in pathDict:

        raise KeyError("Coil height (y) is missing.")

    if pathDict["x"] <= 0.0 or pathDict["y"] <= 0.0:

        raise ValueError("Coil sizes (x/y) must be positive.")

    f = bundle[bundleDict["shape"]](bundleDict, filamentI)

    points = list()

    x = pathDict["x"]
    y = pathDict["y"]
    z = f[1]
    s = [[1, 1], [-1, 1], [-1, -1], [1, -1]]

    for c in range(len(s)):

        p = np.zeros(3)

        p[0] = s[c][0]*x
        p[1] = s[c][1]*y
        p[2] = z

        points.append(p)

    edges = edgeLoopFromPoints(points, edgeStart)

    return points, edges



path = {"loop": pathLoop, "racetrack": pathRaceTrack, "rectangle": pathRectangle}



# +++ Bundles  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #



def bundlePointN(bundleDict):

    return 1

def bundlePointR(bundleDict):

    return 0.0

def bundlePointZ(bundleDict):

    return 0.0

def bundlePointI(bundleDict, I):

    return I

def bundlePoint(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    """

    b = np.zeros(2)

    return b



def bundleCircleN(bundleDict):

    if not "n" in bundleDict:

        raise KeyError("Number of filaments (n) is missing.")

    return bundleDict["n"]

def bundleCircleR(bundleDict):

    if not "r" in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    return 2.0*bundleDict["r"]

def bundleCircleZ(bundleDict):

    if not "r" in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    return 2.0*bundleDict["r"]

def bundleCircleI(bundleDict, I):

    return I/bundleCircleN(bundleDict)

def bundleCircle(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    n : int, Number of filaments
    r : float, Coil bundle radius
    """

    if not "n" in bundleDict:

        raise KeyError("Number of filaments (n) is missing.")

    if not isinstance(bundleDict["n"], int):

        raise KeyError("Number of filaments (n) needs to be of type int.")

    if not bundleDict["n"] > 0:

        raise ValueError("Number of filaments (n) needs to be larger than 0.")

    if not "r" in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    if bundleDict["r"] <= 0.0:

        raise ValueError("Coil bundle radius (r) must be positive.")

    if not i < bundleDict["n"]:

        raise ValueError("Coil filament index (i) out of range (max: n).")

    r = bundleDict["r"]
    phii = 1.0/bundleDict["n"] * 2.0*m.pi

    b = r * np.array([m.cos(i * phii), m.sin(i * phii)])

    return b



def bundleRectangleN(bundleDict, s=4):

    if not "n" in bundleDict:

        raise KeyError("Number of filaments per side (n) is missing.")

    return s*(bundleDict["n"]-1)

def bundleRectangleR(bundleDict):

    if not "r" in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    return bundleDict["r"]

def bundleRectangleZ(bundleDict):

    if not "z" in bundleDict:

        raise KeyError("Coil bundle height (z) is missing.")

    return bundleDict["z"]

def bundleRectangleI(bundleDict, I):

    return I/bundleRectangleN(bundleDict)

def bundleRectangle(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    n : int, Number of filaments per side
    r : float,  Coil size in radial direction
    z : float,  Coil size in axial direction
    """

    if not "n" in bundleDict:

        raise KeyError("Number of filaments per side (n) is missing.")

    if not isinstance(bundleDict["n"], int):

        raise KeyError("Number of filaments per side (n) needs to be of type int.")

    if not bundleDict["n"] > 1:

        raise ValueError("Number of filaments per side (n) needs to be larger than 1.")

    if not "r" in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    if not "z" in bundleDict:

        raise KeyError("Coil bundle height (z) is missing.")

    if bundleDict["r"] <= 0.0 or bundleDict["z"] <= 0.0:

        raise ValueError("Coil bundle sizes (r/z) must be positive.")

    def N(s): return bundleRectangleN(bundleDict, s)
    def N0(s): return bundleRectangleN(bundleDict, s) + 1

    if not i < N(4):

        raise ValueError("Coil filament index (i) out of range (max: 4*(n-1)).")

    b = -0.5 * np.array([bundleDict["r"], bundleDict["z"]])

    bi = -2.0 * b / (bundleDict["n"] - 1)

    if (i < N0(1)):

        b[0] += bi[0] * i
        b[1] += bi[1] * 0.0

    elif (i > N(1)) and (i < N0(2)):

        b[0] += bi[0] * N(1)
        b[1] += bi[1] * (i - N(1))

    elif (i > N(2)) and (i < N0(3)):

        b[0] += bi[0] * (N(3) - i)
        b[1] += bi[1] * N(1)

    elif (i > N(3)) and (i < N(4)):

        b[0] += bi[0] * 0.0
        b[1] += bi[1] * (i - N(3))

    return b



bundleN = {"point": bundlePointN, "circle": bundleCircleN, "rectangle": bundleRectangleN}

bundleI = {"point": bundlePointI, "circle": bundleCircleI, "rectangle": bundleRectangleI}

bundleR = {"point": bundlePointR, "circle": bundleCircleR, "rectangle": bundleRectangleR}

bundleZ = {"point": bundlePointZ, "circle": bundleCircleZ, "rectangle": bundleRectangleZ}

bundle = {"point": bundlePoint, "circle": bundleCircle, "rectangle": bundleRectangle}



# +++ Writing  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

def writeCoilFeatureEdgeMeshes(case, coils):

    for n in coils.keys():

        writeCoilFeatureEdgeMesh(case, coils[n].name,
                                 coils[n].points, coils[n].edges)

def writeCoilFeatureEdgeMesh(case, name, points, edges):

    path = case + "/" + "constant" + "/" + "featureEdgeMesh"
    if not os.path.exists(path): os.makedirs(path)
    fullpath = path + "/" + name + ".eMesh"

    with open(fullpath, "w") as f:

        # Define short indented line with line break
        def ind(iL, cS, eS="\n"): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, "featureEdgeMesh"))

        f.write(ind(0, "// points:\n"))
        f.write(ind(0, str(len(points))))
        f.write(ind(0, "("))

        for p in points:

            f.write(ind(1, "(" + str(p[0]) + " " + str(p[1]) + " " + str(p[2])+")"))

        f.write(ind(0, ")\n"))

        f.write(ind(0, "// edges:\n"))
        f.write(ind(0, str(len(edges))))
        f.write(ind(0, "("))

        for e in edges:

            f.write(ind(1, "(" + str(e[0]) + " "+str(e[1]) + ")"))

        f.write(ind(0, ")\n"))

        f.write(objectFooter())



def writeEdgeBiotSavartProperties(case, coils, nNonOrth=10):

    path = case + "/" + "constant"
    if not os.path.exists(path): os.makedirs(path)
    name = "edgeBiotSavartProperties"
    fullpath = path + "/" + name

    with open(fullpath, "w") as f:

        # Define short indented line with line break
        def ind(iL, cS, eS="\n"): return objectIndent(cS + eS, iLevel=iL)

        # Define write for boolean-strings for OpenFOAM
        def bstr(b): return "true" if b else "false"

        f.write(objectHeader(name, "dictionary"))

        f.write(ind(0, "nNonOrthogonalCorrectors    " + str(nNonOrth) + ";\n"))

        f.write(ind(0, "inductors"))
        f.write(ind(0, "{"))

        for i in range(len(coils)):

            name      = coils[i].name
            reverse   = coils[i].reverse
            current   = coils[i].filamentCurrent
            filaments = coils[i].filaments
            phase     = coils[i].phase

            f.write(ind(1, name))
            f.write(ind(1, "{"))

            f.write(ind(2, "file       " + "\"" + name + ".eMesh\"" + ";"))
            f.write(ind(2, "reverse    " + bstr(reverse) + ";"))
            f.write(ind(2, "current    " + str(current) + ";"))
            f.write(ind(2, "filaments  " + str(filaments) + ";"))
            f.write(ind(2, "phase      " + str(phase) + ";"))

            f.write(ind(1, "}"))

            if not i == len(coils) - 1: f.write(ind(1, ""))


        f.write(ind(0, "}\n"))

        f.write(objectFooter())



def writeFrequency(case, value):

    path = case + "/" + "constant"
    if not os.path.exists(path): os.makedirs(path)
    name = "f0"
    fullpath = path + "/" + name

    with open(fullpath, "w") as f:

        # Define short indented line with line break
        def ind(iL, cS, eS="\n"): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, "uniformDimensionedScalarField"))

        f.write(ind(0, "dimensions    [0 0 -1 0 0 0 0];\n"))
        f.write(ind(0, "value         " + str(value) + ";\n"))

        f.write(objectFooter())

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

class inductorCoil(object):

    def __init__(self, name, bundleDict, pathDict,
                 reverse, current, phase, *args, **kwargs):

        if not isinstance(name, (str, unicode)):

            raise KeyError("Name must be of type string")

        if not isinstance(reverse, bool):

            raise KeyError("Reverse switch must be of type bool.")

        if not isinstance(current, (int, float)):

            raise KeyError("Current must be of type int or float.")

        if not isinstance(phase, (int, float)):

            raise KeyError("Phase must be of type int or float.")

        self.name = name

        self.bundleDict = bundleDict
        self.pathDict = pathDict

        self.reverse = reverse
        self.current = current
        self.filaments = None
        self.filamentCurrent = None
        self.phase = phase

        self.points = list()
        self.edges = list()

        self.compute(*args, **kwargs)

    # ----------------------------------------------------------------------- #

    def compute(self, *args, **kwargs):

        bundleShape = self.bundleDict["shape"]
        pathShape   = self.pathDict["shape"]

        self.filaments = bundleN[bundleShape](self.bundleDict)
        self.filamentCurrent = bundleI[bundleShape](self.bundleDict, self.current)

        for i in range(self.filaments):

            start = len(self.points)

            p, e = path[pathShape](self.pathDict, self.bundleDict, i, start)

            self.points += p
            self.edges += e

        self.transform(*args, **kwargs)

    # ----------------------------------------------------------------------- #

    def transform(self, **kwargs):

        if "function" in kwargs:

            function = kwargs["function"]

            p = self.points
            for i in range(len(p)): p[i] = function(p[i])

        if "translate" in kwargs:

            translate = kwargs["translate"]

            if not isinstance(translate, (list, np.ndarray)):

                raise KeyError("Translation vector must be a list or an array.")

            if isinstance(translate, list): translate = np.array(translate)

            if not len(translate) == 3:

                raise KeyError("Translation vector must have exactly 3 components.")

            p = self.points
            for i in range(len(p)): p[i] += translate

        if "rotate" in kwargs:

            rotate = kwargs["rotate"]

            if not isinstance(rotate, tuple) and len(rotate) == 2:

                raise KeyError("Rotation data must be a tuple containing axis and angle.")

            rotAxis  = rotate[0]
            rotAngle = rotate[1]/180.0 * m.pi

            if not isinstance(rotAxis, (list, np.ndarray)):

                raise KeyError("Rotation vector must be a list or an array.")

            if isinstance(rotAxis,list): rotAxis = np.array(rotAxis)

            if not len(rotAxis) == 3:

                raise KeyError("Rotation vector must have exactly 3 components.")

            rotM = rotationMatrix(rotAxis, rotAngle)

            p = self.points
            for i in range(len(p)): p[i] = np.dot(rotM, p[i])

        if "scale" in kwargs:

            scale = kwargs["scale"]

            if not isinstance(scale, (int, float, list, np.ndarray)):

                raise KeyError("Scalae factor must be of type int, float, list or array.")

            if isinstance(scale, (int, float)): scale = np.array(3 * [scale])

            if isinstance(scale, list): scale = np.array(scale)

            if not len(scale) == 3:

                raise KeyError("Scaling vector must have exactly 3 components.")

            p = self.points
            for i in range(len(p)): p[i] *= scale

    # ----------------------------------------------------------------------- #

    def printData(self):

        print("name:", self.name)
        print("reverse:", reverse)
        print("current:", current)
        print("filaments", filaments)
        print("filamentCurrent", filamentCurrent)
        print("phase:", phase)

        print("pathDict:", self.pathDict)
        print("bundleDict:", self.bundleDict)

        print("points:", self.points)
        print("edges:", self.edges)



class inductorCoils(dict):

    def __setitem__(self, key, item):

        self.__dict__[key] = item

    # ----------------------------------------------------------------------- #

    def __getitem__(self, key):

        return self.__dict__[key]

    # ----------------------------------------------------------------------- #

    def __repr__(self):

        return repr(self.__dict__)

    # ----------------------------------------------------------------------- #

    def __len__(self):

        return len(self.__dict__)

    # ----------------------------------------------------------------------- #

    def __delitem__(self, key):

        del self.__dict__[key]

    # ----------------------------------------------------------------------- #

    def clear(self):

        return self.__dict__.clear()

    # ----------------------------------------------------------------------- #

    def copy(self):

        return self.__dict__.copy()

    # ----------------------------------------------------------------------- #

    def has_key(self, k):

        return self.__dict__.has_key(k)

    # ----------------------------------------------------------------------- #

    def pop(self, k, d=None):

        return self.__dict__.pop(k, d)

    # ----------------------------------------------------------------------- #

    def update(self, *args, **kwargs):

        return self.__dict__.update(*args, **kwargs)

    # ----------------------------------------------------------------------- #

    def keys(self):

        return self.__dict__.keys()

    # ----------------------------------------------------------------------- #

    def values(self):

        return self.__dict__.values()

    # ----------------------------------------------------------------------- #

    def items(self):

        return self.__dict__.items()

    # ----------------------------------------------------------------------- #

    def pop(self, *args):

        return self.__dict__.pop(*args)

    # ----------------------------------------------------------------------- #

    def __cmp__(self, dict):

        return cmp(self.__dict__, dict)

    # ----------------------------------------------------------------------- #

    def __contains__(self, item):

        return item in self.__dict__

    # ----------------------------------------------------------------------- #

    def __iter__(self):

        return iter(self.__dict__)

    # ----------------------------------------------------------------------- #

    def __init__(self, shape, *args, **kwargs):

        super(inductorCoils, self).__init__()

        self._makeCoils[shape](self, *args, **kwargs)

    # ----------------------------------------------------------------------- #

    def _makeCoilsArray(self, name, bundleDict, pathDict,
                        current, n, step, reverse=False,
                        origin=np.zeros(3), axis=2, scale=1.0, **kwargs):

        for i in range(n):

            namei = name + str(i)
            reversei = reverse
            currenti = current
            phasei = 0.0

            translation = np.array(origin).copy()
            translation[axis] += i*step

            self[i] = inductorCoil(namei, bundleDict, pathDict,
                                   reversei, currenti, phasei,
                                   translate=translation, scale=scale)

    def _makeCoilsTMF(self, name, bundleDict, pathDict,
                      current, n, step, reverse=False,
                      origin=np.zeros(3), axis=2, scale=1.0, period=6,
                      stepdouble=None, **kwargs):

        nj = 1

        if stepdouble: nj = 2

        for i in range(n):

            for j in range(nj):

                ij = i + n*j

                namei = name + str(ij)
                reversei = reverse
                currenti = current
                phasei = i*1.0/period * 360.0

                translation = np.array(origin).copy()
                translation[axis] += i*step

                newPathDict = pathDict.copy()

                if j == 1: newPathDict["r"] += stepdouble

                self[ij] = inductorCoil(namei, bundleDict, newPathDict,
                                        reversei, currenti, phasei,
                                        translate=translation, scale=scale)

    def _makeCoilsRMF(self, name, bundleDict, pathDict,
                      current, n, step, reverse=False,
                      origin=np.zeros(3), axis=2, scale=1.0, **kwargs):

        for i in range(n):

            namei = name + str(i)
            reversei = reverse
            currenti = current
            phasei = i*1.0/n * 360.0

            self[i] = inductorCoil(namei, bundleDict, pathDict,
                                   reversei, currenti, phasei)

            plane = [0, 1, 2]; plane.pop(axis)

            rotation = (np.zeros(3), 90.0)
            rotation[0][plane[1]] = 1.0
            self[i].transform(rotate=rotation)

            rotation = (np.zeros(3), 90.0)
            rotation[0][plane[0]] = 1.0
            self[i].transform(rotate=rotation)

            rotation = (np.zeros(3), phasei)
            rotation[0][axis] = 1.0
            self[i].transform(rotate=rotation)

            translation = np.array(origin).copy()
            translation[plane[0]] += step * m.cos(phasei/180.0 * m.pi)
            translation[plane[1]] += step * m.sin(phasei/180.0 * m.pi)
            self[i].transform(translate=translation)

            self[i].transform(scale=scale)

    # ----------------------------------------------------------------------- #

    _makeCoils = {"ARRAY": _makeCoilsArray,
                  "TMF": _makeCoilsTMF,
                  "RMF": _makeCoilsRMF}

    # ----------------------------------------------------------------------- #

    def compute(self, *args, **kwargs):

        for k in self.keys():

            self[k].compute(*args, **kwargs)

    # ----------------------------------------------------------------------- #

    def transform(self, *args, **kwargs):

        for k in self.keys():

            self[k].transform(*args, **kwargs)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

