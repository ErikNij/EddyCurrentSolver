#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module for writing field data for foam from raw data files
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# TODO [High]: Allow other primitive patch types than 'fixedValue'

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

from foamTools.classInfo import (fieldClassToAtomicTypeSize,
                                 patchBaseTypeToPrimitiveValued)
from foamTools.ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fieldDimensions(dim):

    # return : dimensions string
    #
    # dim    : field object dimensions string (e.g. '0 1 -1 0 0 0 0')

    # Define short objectIndented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Assemble dimensions string
    r  = ''
    r += i(0, 'dimensions    [' + str(dim).strip("[]") + ']' + ';')
    r += i(0, '')

    return r



def fieldInternal(dType, dSize, path, offset, comm):

    # return : data string
    #
    # dType  : field data type
    # dSize  : field data type size
    # path   : path to raw internal data
    # offset : column offset for raw data
    # comm   : comments prefix

    # Define short objectIndented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Assemble internal data string
    r  = ''
    r += i(0, 'internalField nonuniform List<' + dType + '>')
    r += i(0, str(sum(1 for l in open(path) if not re.match(comm,l))))
    r += i(0, '(')
    with open(path,'r') as s:
        for l in s:
            if dSize==1:
                r += i(1, ' '.join(l.split()[offset:offset+dSize]))
            if dSize==3:
                r += i(1, '(' + ' '.join(l.split()[offset:offset+dSize]) + ')')
    r += i(0, ')' + ';')
    r += i(0, '')

    return r



def fieldBoundary(dType, dSize, pInfoList):

    # return    : data string
    #
    # dType     : field data type
    # dSize     : field data type size
    # pInfoList : patch info list

    # Define short objectIndented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Include boundary entry
    r  = ''
    r += i(0, 'boundaryField')
    r += i(0, '{')
    r += i(0, '')

    # Assemble boundary data string
    for pInfo in pInfoList:
        p = fieldPatch(dType, dSize, pInfo)
        for l in p.splitlines(True):
            r += i(1,l,eS='')

    r += i(0, '}')

    return r



def fieldPatch(dType, dSize, pInfo):

    # return    : data string
    #
    # dType : field data type
    # dSize : field data type size
    # pInfo : patch info

    # Read patch info
    pName     = pInfo[0]
    pBaseType = pInfo[1]
    pData     = pInfo[2]

    # Get primitive type and valued boolean
    pType, pValued = patchBaseTypeToPrimitiveValued(pBaseType)

    # Define short objectIndented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Include patch name
    r  = ''
    r += i(0, pName)
    r += i(0, '{')

    # Process all patches
    if True:

        # Include patch type
        r += i(1, 'type ' + pType + ';')

    # Process only valued patches
    if pValued:

        # Read patch data info
        path   = pData[0]
        offset = pData[1]
        comm   = pData[2]

        # Assemble patch data string
        r += i(1, 'value nonuniform List<' + dType + '>')
        r += i(1, str(sum(1 for l in open(path) if not re.match(comm,l))))
        r += i(1, '(')
        with open(path,'r') as s:
            for l in s:
                if dSize==1:
                    r += i(2, ' '.join(l.split()[offset:offset+dSize]))
                if dSize==3:
                    r += i(2, '(' + ' '.join(l.split()[offset:offset+dSize]) + ')')
        r += i(1, ')' + ';')

    r += i(0, '}')
    r += i(0, '')

    return r



def fieldWriteFromRaw(time, fInfo, rInfo, pInfoList):

    # time        : time instant (time directory)
    # fInfo       : field information
    #               (0: name, 1: class, 2: dimensions)
    # rInfo       : raw information for internal field data
    #               (0: path, 1: offset, 3: comments)
    # pInfoList   : patch info
    #               (0: name, 1: base type, 2: (raw info))
    #  (raw info) : raw info for patch
    #               (0: path, 1: offset, 3: comments)

    # Read field info
    fName       = fInfo[0]
    fClass      = fInfo[1]
    fDimensions = fInfo[2]

    # Get atomic data type and its size from object class
    fType, fSize = fieldClassToAtomicTypeSize(fClass)

    # Assemble relative path to field file
    fPath = time + '/' + fName

    # Read raw info for internal data
    rawPath    = rInfo[0]
    rawOffset  = rInfo[1]
    rawComment = rInfo[2]

    # Open field file
    with open(fPath,'wx') as f:

        # Define short write
        def fw(cS): f.write(cS)

        # Write header
        fw(objectHeader(fName, fClass, time))

        # Write dimensions
        fw(fieldDimensions(fDimensions))

        # Write internal field
        fw(fieldInternal(fType, fSize, rawPath, rawOffset, rawComment))

        # Write boundary field
        fw(fieldBoundary(fType, fSize, pInfoList))

        # Write footer
        fw(objectFooter())

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

