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

from foamTools.ioInfo import ioBase

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

class setSetBatch(object):

    def __init__(self, filePath=None, maxLength=60):

        self.nameStringMaxLenth = maxLength

        self.io = ioBase(filePath)

        self.io.write()

    # ----------------------------------------------------------------------- #

    def _topoSetString(self, string):

        return string.ljust(8, " ")

    # ----------------------------------------------------------------------- #

    def _nameString(self, string):

        return string.ljust(self.nameStringMaxLenth, " ")

    # ----------------------------------------------------------------------- #

    def _actionString(self, string):

        return string.ljust(6, " ")

    # ----------------------------------------------------------------------- #

    def _writeTopoSetString(self, string):

        self.io.write(self._topoSetString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _writeNameString(self, string):

        self.io.write(self._nameString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _writeActionString(self, string):

        self.io.write(self._actionString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _doAction(self, topo, name, action, source=False, skip=False):

        self._writeTopoSetString(topo + "Set")
        self._writeNameString(topo + "Set" + "_" + name)
        self._writeActionString(action)

        if source:

            self.io.write(source)

        else:

            self.io.write()

        if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doRemove(self, topo, name, skip=False):

        self._doAction(topo, name, "remove", source=False, skip=skip)

    # ----------------------------------------------------------------------- #

    def _doClear(self, topo, name, rmOld=False, skip=True):

        self._doAction(topo, name, "clear", source=False, skip=False)

        if rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doInvert(self, topo, name, rmOld=False, skip=True):

        self._doAction(topo, name, "invert", source=False, skip=False)

        if rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doNew(self, topo, name, source=False, skip=True):

        self._doAction(topo, name, "new", source=source, skip=skip)

    # ----------------------------------------------------------------------- #

    def _doAdd(self, topo, name, source, rmOld=False, skip=True):

        self._doAction(topo, name, "add", source=source, skip=False)

        if rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doDelete(self, topo, name, source, rmOld=False, skip=True):

        self._doAction(topo, name, "delete", source=source, skip=False)

        if rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doAll(self, topo, name, skip=True):

        self._doNew(topo, name, source=False, skip=False)
        self._doInvert(topo, name, rmOld=False, skip=skip)

    # ----------------------------------------------------------------------- #

    def _doInvSetToTopo(self, topo, name, add, rmOld=False, skip=True):

        aList = list(add) if isinstance(add, list) else list([add])

        count = len(aList)

        for a in aList:

            source  = "setTo" + topo.capitalize() + " "
            source += topo + "Set" + "_" + a

            self._doAdd(topo, name, source=source, rmOld=False, skip=False)

        self._doInvert(topo, name, rmOld=False, skip=False)

        if count > 1 and rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doAddSetToTopo(self, topo, name, add, new=True, rmOld=False, skip=True):

        aList = list(add) if isinstance(add, list) else list([add])

        count = len(aList)

        if new:

            self._doNew(topo, name, source=False, skip=False)

        for a in aList:

            source  = "setTo" + topo.capitalize() + " "
            source += topo + "Set" + "_" + a

            self._doAdd(topo, name, source=source, rmOld=False, skip=False)

        if count > 1 and rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doDeleteSetToTopo(self, topo, name, delete, rmOld=False, skip=True):

        dList = list(delete) if isinstance(delete, list) else list([delete])

        count = len(dList)

        for d in dList:

            source  = "setTo" + topo.capitalize() + " "
            source += topo + "Set" + "_" + d

            self._doDelete(topo, name, source=source, rmOld=False, skip=False)

        if count > 1 and rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def _doSetToTopo(self, topo, name, add, delete=None, new=True, invert=False,
                        rmOld=False, skip=True):

        count = len(add) if isinstance(add, list) else 1

        self._doAddSetToTopo(topo, name, add=add, new=new,
                             rmOld=False, skip=False)

        if delete:

            count += len(delete) if isinstance(delete, list) else 1

            self._doDeleteSetToTopo(topo, name, delete=delete,
                                    rmOld=False, skip=False)

        if invert:

            count += 1

            self._doInvert(topo, name, rmOld=False, skip=False)

        if count > 1 and rmOld:

            self._doRemove(topo, name + "_" + "old", skip=skip)

        else:

            if skip: self.io.write()


    # ----------------------------------------------------------------------- #

    def _doSurfaceToTopo(self, topo, name, stl,
                         tol=1e-4, inside=False, outside=False, skip=True):

        # Define write for boolean-strings for OpenFOAM
        def bstr(b): return "true" if b else "false"

        source  = "surfaceToPoint" + " " + "\"" + stl + "\"" + " "
        source += str(tol) + " " + bstr(inside) + " " + bstr(outside)

        self._doNew("point", name, source=source, skip=False)

        if not topo == "point":

            source  = "pointTo" + topo.capitalize() + " "
            source += "pointSet" + "_" + name + " " + "all"

            self._doNew(topo, name, source=source, skip=False)

            self._doRemove("point", name, skip=False)

        if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    def __doCellBoundaryToTopo(self, topo, name, add, delete=None, new=True,
                               rmOld=False, skip=True):

        aList = list(add) if isinstance(add, list) else list([add])

        count = len(aList)

        if topo != "face":

            raise NotImplementedError("Only available for faces.")

        else:

            dList = []

            if delete:

                dList = list(delete) if isinstance(delete, list) else list([delete])

            adList = aList + dList

            for i, ad in enumerate(adList):

                namei = name + str(i)

                self._doNew(topo, namei, source=False, skip=False)

                source  = "cellTo" + topo.capitalize() + " "
                source += "cellSet" + "_" + ad + " " + "all"

                self._doAdd(topo, namei, source=source,
                            rmOld=False, skip=False)

                source  = "cellTo" + topo.capitalize() + " "
                source += "cellSet" + "_" + ad + " " + "both"

                self._doDelete(topo, namei, source=source,
                               rmOld=rmOld, skip=False)

            if skip: self.io.write()

            if new:

                self._doNew(topo, name, source=False, skip=False)

            for i, a in enumerate(aList):

                namei = name + str(i)

                self._doAddSetToTopo(topo, name, add=[namei], new=False,
                                     rmOld=False, skip=False)

            for i, d in enumerate(dList):

                namei = name + str(len(aList) + i)

                self._doDeleteSetToTopo(topo, name, delete=[namei],
                                     rmOld=False, skip=False)

            if skip: self.io.write()

            for i, a in enumerate(adList):

                namei = name + str(i)

                self._doRemove(topo, namei, skip=False)

            if count > 1 and rmOld:

                self._doRemove(topo, name + "_" + "old", skip=skip)

            else:

                if skip: self.io.write()

    # ----------------------------------------------------------------------- #

    _doPoint = {"remove": _doRemove, "clear": _doClear, "invert": _doInvert,
               "new": _doNew, "add": _doAdd, "delete": _doDelete,
               "all": _doAll,
               "invSetToTopo": _doInvSetToTopo,
               "addSetToTopo": _doAddSetToTopo,
               "deleteSetToTopo": _doDeleteSetToTopo,
               "setToTopo": _doSetToTopo,
               "surfaceToTopo": _doSurfaceToTopo,
               "manual": _doAction}

    _doFace = {"remove": _doRemove, "clear": _doClear, "invert": _doInvert,
               "new": _doNew, "add": _doAdd, "delete": _doDelete,
               "all": _doAll,
               "invSetToTopo": _doInvSetToTopo,
               "addSetToTopo": _doAddSetToTopo,
               "deleteSetToTopo": _doDeleteSetToTopo,
               "setToTopo": _doSetToTopo,
               "surfaceToTopo": _doSurfaceToTopo,
               "cellBoundaryToTopo" : __doCellBoundaryToTopo,
               "manual": _doAction}

    _doCell = {"remove": _doRemove, "clear": _doClear, "invert": _doInvert,
               "new": _doNew, "add": _doAdd, "delete": _doDelete,
               "all": _doAll,
               "invSetToTopo": _doInvSetToTopo,
               "addSetToTopo": _doAddSetToTopo,
               "deleteSetToTopo": _doDeleteSetToTopo,
               "setToTopo": _doSetToTopo,
               "surfaceToTopo": _doSurfaceToTopo,
               "manual": _doAction}

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

    def quit(self):

        self.group("End")

        self.io.write("quit")

        self.io.close()

    # ----------------------------------------------------------------------- #

    def write(self, string, end="\n"):

        self.io.write(string, end=end)

    # ----------------------------------------------------------------------- #

    def group(self, string):

        if self.io.filePath and not self.io.fileObject:

            self.io.open()

        self.write("# " + string + "\n")
        return True

    # ----------------------------------------------------------------------- #

    def pointSet(self, name, do, *args, **kwargs):

# TODO: Prevent name from containing spaces

        self._doPoint[do](self, "point", name, *args, **kwargs)

    # ----------------------------------------------------------------------- #

    def faceSet(self, name, do, *args, **kwargs):

# TODO: Prevent name from containing spaces

        self._doFace[do](self, "face", name, *args, **kwargs)

    # ----------------------------------------------------------------------- #

    def cellSet(self, name, do, *args, **kwargs):

# TODO: Prevent name from containing spaces

        self._doCell[do](self, "cell", name, *args, **kwargs)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

