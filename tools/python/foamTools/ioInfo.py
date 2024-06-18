#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# TODO [High]: Allow other primitive patch types than "fixedValue"

# TODO [Low]: Rework function descriptions

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

import re

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def objectIndent(cString, iLevel=0, iChar=" ", iCount=4):

    # return  : indented content string
    #
    # cString : content string to indent
    # iLevel  : indentation level
    # iChar   : indentation character
    # iCount  : indentation character count

    return iLevel * iChar * iCount + cString



def objectHeader(name, cl, time=""):

    # return : header string
    #
    # time   : time instant
    # name   : name of foam object
    # cl     : class of foam object (e.g. volScalarField)

    # Define short indented line with line break
    def i(iL,cS,eS="\n"): return objectIndent(cS + eS,iLevel=iL)

    # Assemble header string
    r  = ""
    r += i(0, "/*--------------------------------*- C++ -*----------------------------------*\\")
    r += i(0, "| =========                 |                                                 |")
    r += i(0, "| \\\\      /  F ield         | foam-extend: Open Source CFD                    |")
    r += i(0, "|  \\\\    /   O peration     | Version:     4.0                                |")
    r += i(0, "|   \\\\  /    A nd           | Web:         http://www.extend-project.de       |")
    r += i(0, "|    \\\\/     M anipulation  | For copyright notice see file Copyright         |")
    r += i(0, "\*---------------------------------------------------------------------------*/")
    r += i(0, "FoamFile")
    r += i(0, "{")
    r += i(1, "version     2.0;")
    r += i(1, "format      ascii;")
    r += i(1, "class       " + cl + ";")
    if time: r += i(1, "location    "" + time + """ + ";")
    r += i(1, "object      " + name + ";")
    r += i(0, "}")
    r += i(0, "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
    r += i(0, "")

    return r



def objectFooter():

    # return : footer string

    # Define short indented line with line break
    def i(iL,cS,eS="\n"): return objectIndent(cS + eS,iLevel=iL)

    # Assemble footer string
    r  = ""
    r += i(0, "// ************************************************************************* //")
    r += i(0, "")

    return r

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

class ioBase(object):

    # ----------------------------------------------------------------------- #

    def __init__(self, filePath=None):

        self.filePath = filePath
        self.fileObject = None

        self.indentLevel = 0

    # ----------------------------------------------------------------------- #

    def _indent(self, level, string, end):

        if not (isinstance(level, int) and isinstance(string, str)
                and isinstance(end, str)):

            raise KeyError("Argument types are wrong.")

        return objectIndent(string + end, iLevel=level)

    # ----------------------------------------------------------------------- #

    def _write(self, string):

        if not isinstance(string, str):

            raise KeyError("Argument types are wrong.")

        if not self.filePath:

            sys.stdout.write(string)

        elif self.fileObject:

            self.fileObject.write(string)

        else:

            raise IOError("File" + " \"" + self.filePath + "\" "+
                          "is not open.")

    # ----------------------------------------------------------------------- #

    def open(self):

        if self.filePath:

            if not self.fileObject:

                self.fileObject = open(os.path.realpath(self.filePath),"w")

                self.fileObject.write("")

            else:

                raise IOError("File" + " \"" + self.filePath + "\" "+
                            "is already open.")

    # ----------------------------------------------------------------------- #

    def close(self):

        if self.filePath:

            if self.fileObject:

                self.fileObject.close()

                self.fileObject = None

            else:

                raise IOError("File" + " \"" + self.filePath + "\" "+
                            "is not open.")

    # ----------------------------------------------------------------------- #

    def rename(self, filePath=None):

        if self.fileObject:

            self.fileObject.close()

            self.fileObject = None

        self.filePath = filePath

    # ----------------------------------------------------------------------- #

    def indent(self, level=0):

        if isinstance(level, int):

            raise KeyError("Indent level must be an integer larger or equal 0.")

        self.indentLevel = level

    # ----------------------------------------------------------------------- #

    def line(self):

        self._write("\n")

    # ----------------------------------------------------------------------- #

    def write(self, string="", ind=True, end="\n"):

        if not (isinstance(string, str) and isinstance(ind, bool)
                and isinstance(end, str)):

            raise KeyError("Argument types are wrong.")

        level = self.indentLevel
        if not ind: level = 0

        wstr = self._indent(level,str(string), end)
        self._write(wstr)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

