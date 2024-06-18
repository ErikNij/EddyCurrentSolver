#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Module with helper functions to write a blockMeshDict with python
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

import warnings

import math as m
import numpy as np
import scipy as sp
import scipy.optimize as spo

from foamTools.ioInfo import ioBase, objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def expansion_n_ds_root(e, ds, n, l=None):

    ls = ds

    if l: ls /= l

    f = pow(e, 1.0/(n - 1))

    return ls*(1.0 - pow(f, n)) - (1.0 - pow(f, 1))



def expansion_n_ds(e, ds, l=None, tol=1e-16, maxiter=100):

    """
    Return distribution size from expansion ratio and start interval

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        e: Expansion ratio
        ds: Start interval
        l: Length
        tol: Target residual for root search
        maxiter: Maxium iteration count for root search

    Returns:
        Distribution size
    """

    ls = ds

    if l: ls /= l

    def N(n): return expansion_n_ds_root(e=e, ds=ls, n=n)

    n0 = 2

    return int(spo.newton(func=N, x0=n0, tol=tol, maxiter=maxiter))



def expansion_f_ds_root(n, ds, f, l=None):

    ls = ds

    if l: ls /= l

    return ls*(1.0 - pow(f, n)) - (1.0 - pow(f, 1))



def expansion_f_ds(n, ds, l=None, tol=1e-16, maxiter=100):

    """
    Return expansion factor from distribution size and start interval

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        ds: Start interval
        l: Length
        tol: Target residual for root search
        maxiter: Maxium iteration count for root search

    Returns:
        Expansion factor
    """

    if isinstance(n, int) and n > 1:

        ls = ds

        if l: ls /= l

        def F(f): return expansion_f_ds_root(n=n, ds=ls, f=f)

        # Make sure we start root search on the right side of the local max
        f0 = 1.0 + 2.0*(pow(1.0/(ls*n), 1.0/(n - 1)) - 1.0)

        return spo.newton(func=F, x0=f0, tol=tol, maxiter=maxiter)

    elif isinstance(n, int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_f_e(n, e):

    """
    Return expansion factor from distribution size and expansion ratio

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        e: Expansion ratio

    Returns:
        Expansion factor
    """

    if isinstance(n, int) and n > 1:

        return pow(e, 1.0/(n - 1))

    elif isinstance(n, int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_e_ds(n, ds, l=None):

    """
    Return expansion ratio from distribution size and start interval

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        ds: Starrt interval

    Returns:
        Expansion ratio
    """

    if isinstance(n, int) and n > 1:

        ls = ds

        if l: ls /= l

        f = expansion_f_ds(n, ls)

        return pow(f, n - 1)

    elif isinstance(n, int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_e_f(n, f):

    """
    Return expansion ratio from distribution size and expansion factor

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        f: Expansion factor

    Returns:
        Expansion ratio
    """

    if isinstance(n, int) and n > 1:

        return pow(f, n - 1)

    elif isinstance(n, int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_ds_e(n, e, l=None):

    """
    Return start interval from distribution size and expansion ratio

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        e: Expansion ratio
        l: Length

    Returns:
        Start interval
    """

    if isinstance(n, int) and n > 1:

        f = expansion_f_e(n, e)

        if f == 1.0:

            ds = 1.0/n

        else:

            ds = (1.0 - pow(f, 1))/(1.0 - pow(f, n))

        if l: ds *= l

        return ds

    else:

        raise ValueError("n must be an integer, which is larger than 1")



def expansion_de_e(n, e, l=None):

    """
    Return end interval from distribution size and expansion ratio

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        e: Expansion ratio
        l: Length

    Returns:
        End interval
    """

    if isinstance(n, int) and n > 1:

        f = expansion_f_e(n, e)

        if f == 1.0:

            de = 1.0/n

        else:

            de = 1.0 - (1.0 - pow(f, n - 1))/(1.0 - pow(f, n))

        if l: de *= l

        return de

    else:

        raise ValueError("n must be an integer, which is larger than 1")



def expansion(n, e, l=None):

    """
    Return all intervals from distribution size and expansion ratio

    Length:               l
    Distribution size:    n
    Expansion ratio:      e  = de/ds
    Expansion factor:     f  = e^(1/n-1)
    Normalized interval:  di = (1 - f^i)/(1 - f^n)
    Norm. first interval: ds = (1 - f^1)/(1 - f^n)
    Norm. last interval:  de = (1 - f^(n-1))/(1 - f^n)
    Interval:             li = di * l

    Args:
        n: Distribution size
        e: Expansion ratio
        l: Length

    Returns:
        All intervals from ds to de
    """

    if isinstance(n, int) and n > 0:

        f = expansion_f_e(n, e)

        t = np.zeros(n+1)

        if f == 1.0:

            for i in range(n+1): t[i] = float(i)/n

        else:

            for i in range(n+1): t[i] = (1.0 - pow(f, i))/(1.0 - pow(f, n))

        if l: t *= l

        return t

    else:

        raise ValueError("n must be an integer, which is larger than 0")

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

class resizeList(list):

    def __init__(self, default=None):

        self.default = default

    # ----------------------------------------------------------------------- #

    def __setitem__(self, key, value):

        if key >= len(self):

            self += [self.default] * (key - len(self) + 1)

        super(resizeList, self).__setitem__(key, value)



class vertices(object):

    # ----------------------------------------------------------------------- #

    def __init__(self, blockMeshDictRef, nmax=1000):

        self.debug = blockMeshDictRef.debug
        self.mesh = blockMeshDictRef.mesh
        self.io = blockMeshDictRef.io

        self.nmax = nmax

        self.labels = list()
        self.labelIndex  = resizeList()

        self.points = list()

    # ----------------------------------------------------------------------- #

    def _meshPoint(self, point, dir=0):

        if len(point) == 3 and self.mesh["dim"] == 3:

            return point

        if len(point) == 2 and self.mesh["dim"] == 2:

            if dir not in [0, 1]: raise ValueError()

            p = [0.0, 0.0, 0.0]

            n = self.mesh["normal"]

            t = [0, 1, 2]; t.pop(n)

            # 2D mesh transformation
            if not self.mesh["wedge"]:

                delta = self.mesh["thickness"] / 2.0

                p[n] = -delta if dir == 0 else delta

                if n == 1: p[n] *= -1.0

                for i in [0, 1]: p[t[i]] = point[i]

                return p

            # 2D axial mesh transformation
            else:

                angle = m.pi/180.0 * self.mesh["angle"] / 2.0

                def rc(r): return r * m.cos(angle)
                def rs(r): return r * m.sin(angle)

                p[n] = -rs(point[0]) if dir == 0 else rs(point[0])

                if n == 1: p[n] *= -1.0

                p[t[0]] = rc(point[0])
                p[t[1]] = point[1]

                return p

        else:

            raise ValueError()

    # ----------------------------------------------------------------------- #

    def set(self, label, point):

        if not (isinstance(label, int) and isinstance(point, (list, np.ndarray))):

            raise KeyError()

        if label > self.nmax:

            raise ValueError("Vertice label "
                + str(label) + ") out of range (max: "
                + str(self.nmax) + ").")

        if not len(point) == self.mesh["dim"]:

            raise ValueError("Vertice (label "
                + str(label) + "), point dimension ("
                + str(len(point)) + ") mismatch (dim: "
                + str(self.mesh["dim"]) + ").")

        if self.mesh["dim"] == 3:

            label = [label]

        else:

            label = [label, self.nmax + 1 + label]

        for i, l in enumerate(label):

            # Vertice already exists
            try:

                vertice = self.labels.index(l)

                self.points[vertice] = self._meshPoint(point, i)

            # New vertice
            except:

                vertice = len(self.labels)
                self.labels.append(l)
                self.labelIndex[l] = vertice

                self.points.append(self._meshPoint(point, i))

    # ----------------------------------------------------------------------- #

    def copyTranslate(self, par1, verticeLabels, direction):

        if not (isinstance(par1, int) and isinstance(direction, list)
                and len(direction) == 3):

            raise KeyError()

        if isinstance(verticeLabels, int) == int:
            verticeLabels = [verticeLabels]

        if not isinstance(verticeLabels, list): raise KeyError()

        for verticeLabel in verticeLabels:

            originPoint = self.points[self.labelIndex[verticeLabel]]

            translatedPoint = \
                [originPoint[i] + direction[i] for i in range(3)]

            if len(verticeLabels) == 1:
                label = par1

            else:
                label = verticeLabel + par1

            self.set(label, translatedPoint)

    # ----------------------------------------------------------------------- #

    def write(self):

        for vertice, verticeLabel in enumerate(self.labels):

            wstr  = "( "
            for component in self.points[vertice]:
                wstr += str(float(component)) + " "
            wstr += ")"

        # ------------------------------------------------------------------- #

            wstr += " // # " + str(verticeLabel)

        # ------------------------------------------------------------------- #

            self.io.write(wstr)



class blocks(object):

    class cTopo(object):

        originVertice = 0

        base = [i for i in range(3)]

        baseVertices = [[0,1],[0,3],[0,4]]

        baseFaces = [[0,1],[2,3],[4,5]]

        baseEdges = [[0,1,2,3],[4,5,6,7],[8,9,10,11]]

        baseTags = {"x": 0, "y": 1, "z": 2}

        baseGetTags = {0: "x", 1: "y", 2: "z"}

        vertices = [i for i in range(8)]

        verticeBaseBase = [[None,0,[0,1],1,2,[0,2],[0,1,2],[1,2]],
                           [0,None,1,[0,1],[0,2],2,[1,2],[0,1,2]],
                           [[0,1],1,None,0,[0,1,2],[1,2],2,[0,2]],
                           [1,[0,1],0,None,[1,2],[0,1,2],[0,2],2],
                           [2,[0,2],[0,1,2],[1,2],None,0,[0,1],1],
                           [[0,2],2,[1,2],[0,1,2],0,None,1,[0,1]],
                           [[0,1,2],[1,2],2,[0,2],[0,1],1,None,0],
                           [[1,2],[0,1,2],[0,2],2,1,[0,1],0,None]]

        verticeBaseSign = [[None,True,[True,True],True,
                            True,[True,True],[True,True,True],[True,True]],
                           [False,None,True,[False,True],
                            [False,True],True,[True,True],[False,True,True]],
                           [[False,False],False,None,False,[False,False,True],
                            [False,True],True,[False,True]],
                           [False,[True,False],True,None,
                            [False,True],[True,False,True],[True,True],True],
                           [False,[True,False],[True,True,False],[True,False],
                            None,True,[True,True],True],
                           [[False,False],False,[True,False],
                            [False,True,False],False,None,True,[False,True]],
                           [[False,False,False],[False,False],False,
                            [False,False],[False,False],False,None,False],
                           [[False,False],[True,False,False],[True,False],
                            False,False,[True,False],True,None]]

        faces = [i for i in range(6)]

        faceBase = [0,0,1,1,2,2]

        faceOpposites = [1,0,3,2,5,4]

        faceVertices = [[0,3,7,4],[1,2,6,5],
                        [0,1,5,4],[2,3,7,6],
                        [0,1,2,3],[4,5,6,7]]

        faceEdges = [[4,11,7,8],[5,10,6,9],
                     [0,9,3,8],[1,10,2,11],
                     [0,5,1,4],[3,6,2,7]]

        faceTags = {"x-": 0, "x+": 1,
                    "y-": 2, "y+": 3,
                    "z-": 4, "z+": 5}

        faceGetTags = {0: "x-", 1: "x+",
                       2: "y-", 3: "y+",
                       4: "z-", 5: "z+"}

        faceBaseTags = {"x": [0,1], "y": [2,3], "z": [4,5]}

        faceGetBaseTags = {0: "x", 1: "x",
                           2: "y", 3: "y",
                           4: "z", 5: "z"}

        edges = [i for i in range(12)]

        edgeBase = [0,0,0,0,1,1,1,1,2,2,2,2]

        edgeVertices = [[0,1],[3,2],[7,6],[4,5],
                        [0,3],[1,2],[5,6],[4,7],
                        [0,4],[1,5],[2,6],[3,7]]

        edgeFaces = [[2,4],[3,4],[3,5],[2,5],
                     [0,4],[1,4],[1,5],[0,5],
                     [0,2],[1,2],[1,3],[0,3]]

        # ------------------------------------------------------------------- #

        def __init__(self, blocksRef):

            self.blocks = blocksRef

    # ----------------------------------------------------------------------- #

    class cSync(object):

        def __init__(self, blocksRef):

            self.blocks = blocksRef

            self.debug = blocksRef.debug

        # ------------------------------------------------------------------- #

        def _setProcessedBlock(self, block):

            self.processedBlocks[block] = True

        # ------------------------------------------------------------------- #

        def _getProcessedBlock(self,block):

            return self.processedBlocks[block]

        # ------------------------------------------------------------------- #

        def reset(self):

            self.processedBlocks = [False for l in self.blocks.labels]

            for block, blockLabel in enumerate(self.blocks.labels):

                self.processedBlocks[block] = False

        # ------------------------------------------------------------------- #

        def sync(self, block, blockSyncData, blockSync, base):

            if self.debug:

                debugMsg = "[blocks::_setNeighbourData]"

            faces = self.blocks.topo.faces
            baseFaces = self.blocks.topo.baseFaces[base]
            crossBaseFaces = \
                [f for i, f in enumerate(faces) if i not in baseFaces]

            if self.debug:

                print(debugMsg, "block =", block)
                print(debugMsg, "blockSyncData =", blockSyncData)
                print(debugMsg, "base =", base)
                print(debugMsg, "baseFaces =", baseFaces)
                print(debugMsg, "crossBaseFaces =", crossBaseFaces)
                print()


            # Write data for distribution and simpleGrading data
            if (len(blockSync) == 3):

                blockSyncData[block][base] = blockSync[base]

            # Write data of edgeGrading data
            elif (len(blockSync) == 12):

                pass
# TODO [High]: Generalize WRITE to be used for edgeGrading
# TODO [High]: What if syncing edgeGrading but blockSyncData is simpleGrading
# TODO [High]: Now we need to consider the base direction changes!!!

            self._setProcessedBlock(block)

            # Sync
            for face in crossBaseFaces:

                # Get face vertices
                faceVertices = self.blocks.topo.faceVertices[face]
                faceVerticesBaseBase, faceVerticesBaseSign = \
                    self.blocks._getPathVerticeBase(faceVertices)
                faceVerticeLabels = \
                    [ self.blocks.blockVertices[block][v] \
                        for v in faceVertices]

                # Get next neighbour and exit if there is none
                nextBlock = self.blocks.faceNeighbours[block][face]

                if self.debug:

                    print(debugMsg, "face =", face)
                    print(debugMsg, "faceVertices =", faceVertices)
                    print(debugMsg, "faceVerticesBaseBase =",
                          faceVerticesBaseBase)
                    print(debugMsg, "faceVerticesBaseSign =",
                          faceVerticesBaseSign)
                    print(debugMsg, "faceVerticeLabels =", faceVerticeLabels)
                    print(debugMsg, "nextBlock =", nextBlock)
                    if not nextBlock == None:
                        print(debugMsg, "getProcessedBlock(nextBlock) =",
                              self._getProcessedBlock(nextBlock))
                    print()

                if nextBlock == None or \
                    self._getProcessedBlock(nextBlock) == True:

                    continue

                # Get next face index
                nextFace = self.blocks.faceNeighbours[nextBlock].index(block)

                # Get next block vertices
                nextBlockVertices = self.blocks.blockVertices[nextBlock]

                # Get next face vertices
                nextFaceVertices = \
                    self.blocks.topo.faceVertices[nextFace]
                nextFaceVerticesShifted = \
                    [ nextBlockVertices.index(l) for l in faceVerticeLabels ]
                nextFaceVerticeLabels = \
                    [ nextBlockVertices[v] for v in nextFaceVertices]
                nextFaceVerticesBaseBase, nextFaceVerticesBaseSign = \
                    self.blocks._getPathVerticeBase(nextFaceVerticesShifted)

                # Calculate base transformation map
                nextBaseOrientation = [None for i in self.blocks.topo.base]

                nextBaseOrientation[self.blocks.topo.faceBase[face]] = \
                    self.blocks.topo.faceBase[nextFace]
                nextBaseOrientation[faceVerticesBaseBase[0]] = \
                    nextFaceVerticesBaseBase[0]
                nextBaseOrientation[faceVerticesBaseBase[1]] = \
                    nextFaceVerticesBaseBase[1]

                # Calculate reverse base transformation map
                nextbaseOrientationReverse = \
                    [ nextBaseOrientation.index(b) \
                    for b in self.blocks.topo.base]

                # Assemble block sync list for next block
                # for distribution and simpleGrading data
                if (len(blockSync) == 3):

                    nextBlockSync = \
                        [ blockSync[b] for b in nextbaseOrientationReverse ]

                # Assemble block sync list for next block
                # for edgeGrading data
                elif (len(blockSync) == 12):

                    pass
# TODO [High]: Generalize WRITE to be used for edgeGrading
# TODO [High]: What if syncing edgeGrading but blockSyncData is simpleGrading
# TODO [High]: Now we need to consider the base direction changes!!!

                # Extract next base
                nextBase = nextBaseOrientation[base]

                if self.debug:

                    print(debugMsg, "nextFace =", nextFace)
                    print(debugMsg, "nextFaceVertices =", nextFaceVertices)
                    print(debugMsg, "nextFaceVerticesTransformed =",
                        nextFaceVerticesShifted)
                    print(debugMsg, "nextFaceVerticeLabels =",
                        nextFaceVerticeLabels)
                    print(debugMsg, "nextFaceVerticesBaseBase =",
                        nextFaceVerticesBaseBase)
                    print(debugMsg, "nextFaceVerticesBaseSign =",
                        nextFaceVerticesBaseSign)
                    print(debugMsg, "nextBaseOrientation =",
                        nextBaseOrientation)
                    print(debugMsg, "nextbaseOrientationReverse =",
                        nextbaseOrientationReverse)
                    print(debugMsg, "nextBlockSync =", nextBlockSync)
                    print(debugMsg, "nextBase =", nextBase)
                    print()

                # Recursion
                self.sync(nextBlock, blockSyncData, nextBlockSync, nextBase)

    # ----------------------------------------------------------------------- #

    class cDistribution(cSync):

        def __init__(self, blocksRef):

            super(blocksRef.cDistribution, self).__init__(blocksRef)

        # ------------------------------------------------------------------- #

        def set(self, blockLabels, par1, par2=None):

            if isinstance(blockLabels, int): blockLabels = [ blockLabels ]

            if not isinstance(blockLabels, list): raise KeyError()

            blocks = [ self.blocks.labelIndex[l] for l in blockLabels ]

            for block in blocks:

                # Distribution for all components
                if not par2:

                    if isinstance(par1, list) and len(par1) == 3:

                        blockDistribution = par1

                        # Set distribution for current block
                        self.blocks.distributions[block] = blockDistribution

                        # Sync blocks for all components
                        for base in self.blocks.topo.base:

                            self.reset()
                            self.sync(
                                block,
                                self.blocks.distributions,
                                self.blocks.distributions[block],
                                base
                            )

                    else: raise KeyError()

                # Direction and distribution component
                else:

                    # Direction as base
                    if isinstance(par1, int):

                        base = par1

                    # Direction as baseTag
                    elif isinstance(par1, str):

                        baseTag = par1

                        base = self.blocks.topo.baseTags[baseTag]

                    else: raise KeyError()

                    # Distribution component
                    if isinstance(par2, int):

                        blockDistributionComponent = par2

                    else: raise KeyError()

                    # Set distribution component for current block
                    self.blocks.distributions[block][base] = \
                        blockDistributionComponent

                    # Sync blocks for component
                    self.reset()
                    self.sync(
                        block,
                        self.blocks.distributions,
                        self.blocks.distributions[block],
                        base
                    )

    # ----------------------------------------------------------------------- #

    class cGrading(cSync):

        def __init__(self, blocksRef):

            super(blocksRef.cGrading, self).__init__(blocksRef)

        # ------------------------------------------------------------------- #

        def set(self, blockLabels, par1, par2=None):

            if isinstance(blockLabels, int): blockLabels = [ blockLabels ]

            if not isinstance(blockLabels, list): raise KeyError()

            blocks = [ self.blocks.labelIndex[l] for l in blockLabels ]

            for block in blocks:

                # Grading for all components
                if not par2:

                    # simpleGrading
                    if isinstance(par1, list) and len(par1) == 3:

                        blockSimpleGrading = par1

                        # Set grading for current block
                        self.blocks.gradings[block] = blockSimpleGrading

# TODO [High]: Implement edgeGrading
# TODO [High]: What if syncing edgeGrading but blockSyncData is simpleGrading
# TODO [High]: Now we need to consider the base direction changes!!!

                    else: raise KeyError()

                    # Sync blocks for all components
                    for base in self.blocks.topo.base:

                        self.reset()
                        self.sync(
                            block,
                            self.blocks.gradings,
                            self.blocks.gradings[block],
                            base
                        )

                # Direction and grading component
                else:

                    # Direction as base
                    if isinstance(par1, int):

                        base = par1

                    # Direction as baseTag
                    elif isinstance(par1, str):

                        baseTag = par1

                        base = self.blocks.topo.baseTags[baseTag]

                    else: raise KeyError()

                    # Grading component
                    if True:

                        # simpleGrading
                        if isinstance(par2, float):

                            blockSimpleGradingComponent = par2

                            # Set grading component for current block
                            self.blocks.gradings[block][base] = \
                                blockSimpleGradingComponent

# TODO [High]: Implement edgeGrading
# TODO [High]: What if syncing edgeGrading but blockSyncData is simpleGrading
# TODO [High]: Now we need to consider the base direction changes!!!

                        else: raise KeyError()

                    # Sync blocks in given direction
                    self.reset()
                    self.sync(
                        block,
                        self.blocks.gradings,
                        self.blocks.gradings[block],
                        base
                    )

    # ----------------------------------------------------------------------- #

    # ----------------------------------------------------------------------- #

    def __init__(self, blockMeshDictRef, nmax=1000):

        self.debug = blockMeshDictRef.debug
        self.mesh = blockMeshDictRef.mesh
        self.io = blockMeshDictRef.io
        self.vertices = blockMeshDictRef.vertices

        self.nmax = nmax

        self.topo = self.cTopo(self)
        self.distribution = self.cDistribution(self)
        self.grading = self.cGrading(self)

        self.labels = list()
        self.labelIndex = resizeList()

        self.blockVertices = list()
        self.blockVerticeLabels = list()
        self.distributions = list()
        self.gradings = list()
        self.zones = list()

        self.neighbours = list()
        self.faceNeighbours = list()
        self.edgeNeighbours = list()

    # ----------------------------------------------------------------------- #

    def _getPathVerticeBase(self,vertices):

        bases = []
        signs = []

        for i in range(len(vertices)):

            if i < (len(vertices)-1):

                bases.append(
                    self.topo.verticeBaseBase[vertices[i]][vertices[i+1]]
                )
                signs.append(
                    self.topo.verticeBaseSign[vertices[i]][vertices[i+1]]
                )

            else:

                bases.append(
                    self.topo.verticeBaseBase[
                        vertices[len(vertices)-1]
                    ][
                        vertices[0]
                    ]
                )
                signs.append(
                    self.topo.verticeBaseSign[
                        vertices[len(vertices)-1]
                    ][
                        vertices[0]
                    ]
                )

        return bases, signs


    # ----------------------------------------------------------------------- #

    def _getSharedVertices(self, vertices, otherVertices):

        if not (isinstance(vertices, list) and isinstance(otherVertices, list)):

            raise KeyError()

        return list(set(vertices) & set(otherVertices))

    # ----------------------------------------------------------------------- #

    def _getNeighbours(self, block, otherBlock):

        if self.debug:

            debugMsg = "[blocks::_getNeighbours]"
            print(debugMsg, "block =", block)
            print(debugMsg, "otherBlock =", otherBlock)
            print()

        if not (isinstance(block, int) and isinstance(otherBlock, int)):

            raise KeyError()

        # Get both sets of block vertices
        blockVertices = self.blockVertices[block]
        otherBlockVertices = self.blockVertices[otherBlock]

        # Find shared vertices
        neighbourVertices = \
            self._getSharedVertices(blockVertices, otherBlockVertices)

        if self.debug:

            print(debugMsg, "neighbourVertices =", neighbourVertices)

        # Face neighbours
        if len(neighbourVertices) == 4:

            return True, "face", neighbourVertices

        # Wedge axis face neighbours
        if len(neighbourVertices) == 3:

            return True, "axis", neighbourVertices

        # Edge neighbours
        if len(neighbourVertices) == 2:

            return True, "edge", neighbourVertices

        # No neighbours
        else:

            return False, None, neighbourVertices

    # ----------------------------------------------------------------------- #

    def _setFaceNeighbourData(self, block, otherBlock, neighbourVertices):

        blockVertices = self.blockVertices[block]
        otherBlockVertices = self.blockVertices[otherBlock]

        for face, faceVertices in enumerate(self.topo.faceVertices):

            faceVerticesLabels = \
                [ blockVertices[i] for i in faceVertices ]
            otherFaceVerticesLabels = \
                [ otherBlockVertices[i] for i in faceVertices ]

            sharedFaceVertices = self._getSharedVertices(
                neighbourVertices,
                faceVerticesLabels
            )
            otherSharedFaceVertices = self._getSharedVertices(
                neighbourVertices,
                otherFaceVerticesLabels
            )

            if len(sharedFaceVertices) == 4 \
                or len(sharedFaceVertices) == 3:
                self.faceNeighbours[block][face] = otherBlock
            if len(otherSharedFaceVertices) == 4 \
                or len(otherSharedFaceVertices) == 3:
                self.faceNeighbours[otherBlock][face] = block

    # ----------------------------------------------------------------------- #

    def _setEdgeNeighbourData(self, block, otherBlock, neighbourVertices):

        blockVertices = self.blockVertices[block]
        otherBlockVertices = self.blockVertices[otherBlock]

        for edge, edgeVertices in enumerate(self.topo.edgeVertices):

            edgeVerticesLabels = \
                [ blockVertices[i] for i in edgeVertices ]
            otherEdgeVerticesLabels = \
                [ otherBlockVertices[i] for i in edgeVertices ]

            sharedEdgeVertices = self._getSharedVertices(
                neighbourVertices,
                edgeVerticesLabels
            )
            otherSharedEdgeVertices = self._getSharedVertices(
                neighbourVertices,
                otherEdgeVerticesLabels
            )

            if len(sharedEdgeVertices) == 2:
                try:
                    self.edgeNeighbours[block][edge].index(otherBlock)
                except:
                    if self.edgeNeighbours[block][edge] == None:
                        self.edgeNeighbours[block][edge] = []
                    self.edgeNeighbours[block][edge].append(otherBlock)
            if len(otherSharedEdgeVertices) == 2:
                try:
                    self.edgeNeighbours[otherBlock][edge].index(block)
                except:
                    if self.edgeNeighbours[otherBlock][edge] == None:
                        self.edgeNeighbours[otherBlock][edge] = []
                    self.edgeNeighbours[otherBlock][edge].append(block)

    # ----------------------------------------------------------------------- #

    def _setNeighbourData(self, block):

        if self.debug:

            debugMsg = "[blocks::_setNeighbourData]"

            print(debugMsg, "block =", block)
            print()

        for otherBlock, otherBlockLabel in enumerate(self.labels):

            if not block == otherBlock:

                neighbours, neighbourType, neighbourVertices = \
                    self._getNeighbours(block, otherBlock)

                if self.debug:

                    print(debugMsg, "otherBlock =", otherBlock)
                    print()

                    print(debugMsg, "neighbours =", neighbours)
                    print(debugMsg, "neighbourType =", neighbourType)
                    print(debugMsg, "neighbourVertices =", neighbourVertices)
                    print()

                if neighbourType == "face" \
                    or neighbourType == "axis":

                    try:
                        self.neighbours[block].index(otherBlock)
                    except:
                        if self.neighbours[block] == None:
                            self.neighbours[block] = []
                        self.neighbours[block].append(otherBlock)
                    try:
                        self.neighbours[otherBlock].index(block)
                    except:
                        if self.neighbours[otherBlock] == None:
                            self.neighbours[otherBlock] = []
                        self.neighbours[otherBlock].append(block)

                    # Set all face neighbours
                    self._setFaceNeighbourData(
                        block,
                        otherBlock,
                        neighbourVertices
                    )

                    ## Set all edge neighbours
                    self._setEdgeNeighbourData(
                        block,
                        otherBlock,
                        neighbourVertices
                    )

    # ----------------------------------------------------------------------- #

    def set(
        self, label, verticeLabels,
        distribution=None, grading=None, zone=None
    ):

        if not (isinstance(label, int) and isinstance(verticeLabels, list)
                and (distribution == None or isinstance(distribution, list))
                and (grading == None or (isinstance(grading, list)
                                         and (len(grading) == 3
                                           or len(grading) == 12)))
                and (zone == None or isinstance(zone, str))):

            raise KeyError()

        if label > self.nmax:

            raise ValueError("Block label "
                + str(label) + ") out of range (max: "
                + str(self.nmax) + ").")

        if not len(verticeLabels) == 2**self.mesh["dim"]:

            raise ValueError("Block (label "
                + str(label) + "), vertice label dimension ("
                + str(len(verticeLabels)) + ") mismatch (2^dim: "
                + str(2**self.mesh["dim"]) + ").")

        # Add 2D vertices
        if self.mesh["dim"] == 2:

            verticeLabels += [self.vertices.nmax + 1 + l for l in verticeLabels]

        blockVertices = [self.vertices.labelIndex[l] for l in verticeLabels]

        if distribution:

## TODO: Check and processdistribution  for 2D-meshes
            pass

        else:

            distribution = [1,1,1]

        if grading:

## TODO: Check and processdistribution  for 2D-meshes
            pass

        else:

            grading = [1.0,1.0,1.0]

        # ------------------------------------------------------------------- #

        # Block already exists
        try:

            block = self.labels.index(label)

            self.blockVertices[block] = blockVertices
            self.blockVerticeLabels[block] = verticeLabels
            self.distributions[block] = distribution
            self.gradings[block] = grading
            self.zones[block] = zone
            self.neighbours[block] = None
            self.faceNeighbours[block] = [None for i in range(6)]
            self.edgeNeighbours[block] = [None for i in range(12)]

        # New block
        except:

            block = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = block

            self.blockVertices.append(blockVertices)
            self.blockVerticeLabels.append(verticeLabels)
            self.distributions.append(distribution)
            self.gradings.append(grading)
            self.zones.append(zone)
            self.neighbours.append(None)
            self.faceNeighbours.append([None for i in range(6)])
            self.edgeNeighbours.append([None for i in range(12)])

        self._setNeighbourData(block)

    # ----------------------------------------------------------------------- #

    def copyShiftVerticeLabels(self, shift, blockLabels, verticeLebelShift):

        if not (isinstance(shift, int)
                and (isinstance(verticeLebelShift, int)
                     or isinstance(verticeLebelShift, list))):

            raise KeyError()

# TODO: Is this really necessary?
        #if verticeLebelShift < 2*self.nmax:

            #raise ValueError("Vertice label shift ("
                #+ str(verticeLebelShift) + ") is too small (min: "
                #+ str(2*self.nmax) + ").")

        if isinstance(blockLabels, int):
            blockLabels = [ blockLabels ]

        if not isinstance(blockLabels, list): raise KeyError()

        if isinstance(verticeLebelShift, int):
            verticeLebelShift = [ verticeLebelShift for l in range(8)]

        if not isinstance(verticeLebelShift, list): raise KeyError()

        for blockLabel in blockLabels:

            label = blockLabel + shift

            originVerticeLables = \
                self.blockVerticeLabels[self.labelIndex[blockLabel]]

            shiftedVerticeLables = list(originVerticeLables)

            for l in range(8):

                shiftedVerticeLables[l] = \
                    originVerticeLables[l] + verticeLebelShift[l]

# TODO: How to COPY those values? The below are references...
            originDistribution = self.distributions[blockLabel]
            originGrading = self.gradings[blockLabel]
            originZone = self.zones[blockLabel]
# TODO: How to COPY those values? The above are references...

            self.set(label, shiftedVerticeLables)

    # ----------------------------------------------------------------------- #

    def write(self):

        for block, blockLabel in enumerate(self.labels):

            self.io.write("hex", end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for vertice in self.blockVertices[block]:
                wstr += str(vertice) + " "
            wstr += ")"

            self.io.write(wstr, ind=False, end=" ")

        # ------------------------------------------------------------------- #

            if self.zones[block]:

                self.io.write(str(self.zones[block]), ind=False, end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for distribution in self.distributions[block]:
                wstr += str(distribution) + " "
            wstr += ")"

            self.io.write(wstr, ind=False, end=" ")

        # ------------------------------------------------------------------- #

            if len(self.gradings[block]) == 3:

                self.io.write("simpleGrading", ind=False, end=" ")

            elif len(self.gradings[block]) == 12:

                self.io.write("edgeGrading", ind=False, end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for grading in self.gradings[block]:
                wstr += str(grading) + " "
            wstr += ")"

        # ------------------------------------------------------------------- #

            wstr += " // # " + str(blockLabel)

            self.io.write(wstr, ind=False)

    # ----------------------------------------------------------------------- #

    def printData(self, blockLabels=None):

        # Only print(data for given block
        if blockLabels:

            if isinstance(blockLabels, int):
                blockLabels = [ blockLabels ]

            if not isinstance(blockLabels, list): raise KeyError()

            printBlocks = [ self.labelIndex[l] for l in blockLabels ]

        # Print all blocks
        else:

            printBlocks = [ self.labelIndex[l] for l in self.labels ]

        # ------------------------------------------------------------------- #

        for block in printBlocks:

            print("index:", block)
            print("label:", self.labels[block])
            print("vertices: ", self.blockVertices[block])
            print("verticeLabels: ", self.blockVerticeLabels[block])
            print("distributions: ", self.distributions[block])
            print("gradings: ", self.gradings[block])
            print("zone: ", self.zones[block])
            print("neighbours: ", self.neighbours[block])
            print("faceNeighbours: ", self.faceNeighbours[block])
            print("edgeNeighbours: ", self.edgeNeighbours[block])
            print()



class boundaryFaces(object):

    def __init__(self, blockMeshDictRef):

        self.debug = blockMeshDictRef.debug
        self.mesh = blockMeshDictRef.mesh
        self.io = blockMeshDictRef.io
        self.vertices = blockMeshDictRef.vertices
        self.blocks = blockMeshDictRef.blocks

        self.labels = list()
        self.labelIndex = resizeList()

        self.faces = list()
        self.boundary = list()

    # ----------------------------------------------------------------------- #

    def set(self, label, boundary, par1, par2=None):

        if not isinstance(boundary, str): raise KeyError()

        blockLabels = par1

        if isinstance(blockLabels, int):
            blockLabels = [ blockLabels ]

        if not isinstance(blockLabels, list): raise KeyError()

        # Single
        if not par2 == None:
            faceTag = par2

            face = self.blocks.topo.faceTags[faceTag]

            faceVerticeLabelsList = []

            for blockLabel in blockLabels:

                block = self.blocks.labelIndex[blockLabel]
                blockVerticeLabels = self.blocks.blockVertices[block]

                faceVertices = self.blocks.topo.faceVertices[face]
                faceVerticeLabels = \
                    [ blockVerticeLabels[i] for i in faceVertices ]

                faceVerticeLabelsList.append(faceVerticeLabels)

        elif par2 == None and len(par1) == 4:

            faceVerticeLabels = par1

            faceVerticeLabelsList = [ faceVerticeLabels ]

        else: raise KeyError()

        # ------------------------------------------------------------------- #

        # Boundary face already exists
        try:

            face = self.labels.index(label)

            self.boundary[face] = boundary
            self.faces[face] = faceVerticeLabelsList

        # New boundary face
        except:

            face = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = face

            self.boundary.append(boundary)
            self.faces.append(faceVerticeLabelsList)

    # ----------------------------------------------------------------------- #

    def write(self):

        self.io.write("faces")
        self.io.write("(")
        self.io.indentLevel += 1

        # ------------------------------------------------------------------- #

        for boundaryFace, boundaryFaceLabel in enumerate(self.labels):

            if self.boundary[boundaryFace] == self.io.subdict.boundary:

                for face in self.faces[boundaryFace]:

                    wstr  = "( "
                    for vertice in face:
                        wstr += str(vertice) + " "
                    wstr += ")"
                    self.io.write(wstr)

        # ------------------------------------------------------------------- #

        self.io.indentLevel -= 1
        self.io.write(");")

    # ----------------------------------------------------------------------- #

    def printData(self, boundaryFaceLabels=None):

        # Only print data for given boundary face
        if boundaryFaceLabels:

            if isinstance(boundaryFaceLabels, int):
                boundaryFaceLabels = [ boundaryFaceLabels ]

            if not isinstance(boundaryFaceLabels, list): raise KeyError()

            printBoundaryFaces = \
                [ self.labelIndex[l] for l in boundaryFaceLabels ]

        # Print all boundary faces
        else:

            printBoundaryFaces = [ self.labelIndex[l] for l in self.labels ]

        # ------------------------------------------------------------------- #

        for boundaryFace in printBoundaryFaces:

            print("index:", boundaryFace)
            print("label:", self.labels[boundaryFace])
            print("boundary: ", self.boundary[boundaryFace])
            print("faces: ", self.faces[boundaryFace])
            print()



class blockMeshDict(object):

    class cIo(ioBase):

        class subDictData(object):

            opened = False
            boundary = None

        # ------------------------------------------------------------------- #

        def __init__(self, blockMeshDictRef, filePath=None):

            super(blockMeshDictRef.cIo, self).__init__(filePath)

            self.blockMeshDict = blockMeshDictRef

            self.subdict = self.subDictData()

        # ------------------------------------------------------------------- #

        def boundarySubDictHeader(self, name, typename="empty"):

            self.subdict.boundary = name
            self.write(str(name))
            self.write("{")
            self.indentLevel += 1
            self.write("type" + " " + str(typename) + ";")

        # ------------------------------------------------------------------- #

        def boundarySubDictFooter(self, line = True):

            self.indentLevel -= 1
            self.write("}")
            self.subdict.boundary = None
            if line: self.line()

        # ------------------------------------------------------------------- #

        def boundarySubDict(self, name, typename="empty"):

            if self.subdict.boundary: self.boundarySubDictFooter()
            self.boundarySubDictHeader(name, typename)
            return True

        # ------------------------------------------------------------------- #

        def subDictHeader(self, name):

            self.subdict.opened = True
            self.write(str(name))
            self.write("(")
            self.indentLevel += 1

        # ------------------------------------------------------------------- #

        def subDictFooter(self, line = True):

            self.indentLevel -= 1
            self.write(");")
            self.subdict.opened = False
            if line: self.line()

        # ------------------------------------------------------------------- #

        def subDict(self, name):

            if self.subdict.boundary: self.boundarySubDictFooter(line=False)
            if self.subdict.opened: self.subDictFooter()
            self.subDictHeader(name)
            return True

        # ------------------------------------------------------------------- #

        def header(self, scale):

            wstr = objectHeader("blockMeshDict", "dictionary")
            self.write(wstr)
            wstr = "convertToMeters" + " " + str(scale) + ";"
            self.write(wstr)
            self.line()

        # ------------------------------------------------------------------- #

        def footer(self):

            if self.subdict.opened: self.subDictFooter()
            wstr = objectFooter()
            self.write(wstr)

    # ----------------------------------------------------------------------- #

    def __init__(self, filePath=None,
                 nmaxverts=1000, nmaxblocks=1000,
                 mesh=None, debug=False):

        self.debug = debug

        self._initMesh(mesh)
        self.io = self.cIo(self, filePath)
        self.vertices = vertices(self, nmaxverts)
        self.blocks = blocks(self, nmaxblocks)
        self.boundaryFaces = boundaryFaces(self)

    # ----------------------------------------------------------------------- #

    def _initMesh(self, mesh):

        if not mesh: mesh = dict()

        if not "normal" in mesh: mesh["normal"] = -1

        if not mesh["normal"] in [-1, 0, 1, 2]:

            raise ValueError("Mesh normal must be one of [-1, 0, 1, 2].")

        else:

            if mesh["normal"] == -1:

                if not "dim" in mesh: mesh["dim"] = 3

            else:

                if not "dim" in mesh: mesh["dim"] = 2

        if not mesh["dim"] in [2, 3]:

            raise ValueError("Mesh dimensions must be one of [2, 3].")

        if mesh["dim"] == 3:

            mesh.pop("wedge", None)
            mesh.pop("extent", None)
            mesh.pop("angle", None)
            mesh.pop("thickness", None)

        else:

            if not "wedge" in mesh: mesh["wedge"] = False

            if mesh["wedge"]:

                mesh.pop("thickness", None)

                if "extent" in mesh:

                    mesh["angle"] = mesh["extent"]

                else:

                    if not "angle" in mesh: mesh["angle"] = 5.0

                if mesh["angle"] < 2.0 or mesh["angle"] > 8.0:

                    warnings.warn("Wedge mesh angle should preferably be 5°")

                    if mesh["angle"] <= 0.0:

                        raise ValueError("Wedge mesh angle may not be zero or negative.")

            else:

                mesh.pop("angle", None)

                if "extent" in mesh:

                    mesh["thickness"] = mesh["extent"]

                else:

                    if not "thickness" in mesh: mesh["thickness"] = 1.0

                if mesh["thickness"] <= 0.0:

                    raise ValueError("Mesh thickness may not be zero or negative.")

        self.mesh = mesh

    # ----------------------------------------------------------------------- #

    def rename(self, filePath):

        self.io.rename(filePath)

    # ----------------------------------------------------------------------- #

    def open(self):

        self.io.open()

    # ----------------------------------------------------------------------- #

    def close(self):

        self.io.close()

    # ----------------------------------------------------------------------- #

    def write(self, string, end="\n"):

        self.io.write(string, end=end)

    # ----------------------------------------------------------------------- #

    def header(self, scale):

        self.io.open()

        self.io.header(scale)

    # ----------------------------------------------------------------------- #

    def footer(self):

        self.io.footer()

        self.io.close()

    # ----------------------------------------------------------------------- #

    def subDict(self, name):

        self.io.subDict(name)
        return True

    # ----------------------------------------------------------------------- #

    def boundarySubDict(self, name, typename):

        self.io.boundarySubDict(name, typename)
        return True

    # ----------------------------------------------------------------------- #

    def arc(self, label1, label2, label3, end="\n"):

        vertice1 = self.vertices.labelIndex[label1]
        vertice2 = self.vertices.labelIndex[label2]
        vertice3 = self.vertices.labelIndex[label3]

        wstr  = "arc " + str(vertice1) + " " + str(vertice2)

        wstr += " ( "
        for component in self.vertices.points[vertice3]:
            wstr += str(float(component)) + " "
        wstr += ")"

        self.io.write(wstr, end=end)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

