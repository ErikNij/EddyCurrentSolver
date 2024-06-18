#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# TODO [High]: Allow other primitive patch types than "fixedValue"

# TODO [Low]: Rework function descriptions

# TODO [Low]: Reorganize some functions in modules

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

# pyFoam modules
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import (Vector, Field, Dimension, DictProxy,
                                          TupleProxy, Tensor, SymmTensor,
                                          Unparsed, UnparsedList, Codestream,
                                          DictRedirection, BinaryBlob,
                                          BinaryList, BoolProxy)

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def tryKey(d, k):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : key of dictionary to try

    # TODO [Low]: Improve descriptions

    try:
        assert d[k]==d[k]

    except:
        return False

    return True



def tryMandatory(d, k):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : mandatory key of dictionary

    # TODO [Low]: Improve descriptions

    if tryKey(d, k):
        return d[k]

    else:
        raise KeyError("mandatory key " + k + " is missing")



def tryOrDefault(d, k, vDef):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : key of dictionary to try
    # vDef   : default value of key

    # TODO [Low]: Improve descriptions

    if tryKey(d, k):
        return d[k]

    else:
        return vDef



def readStr(v):

    # return : stripped dictionary value (without "\"")
    #
    # v      : dictionary key value of type str

    # TODO [Low]: Improve descriptions

    if isinstance(v, (str, unicode)):
        return v.strip("\"")
    else:
        try:
            return str(v)
        except:
            raise TypeError("key value " + str(v) + " is no string")



def readInt(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type int

    # TODO [Low]: Improve descriptions

    if isinstance(v, int):
        return v

    if isinstance(v, (str, unicode)):
        return int(readStr(v))

    else:
      raise TypeError("key value " + str(v) + " is not an integer")



def readFloat(v):

    # return : float dictionary value
    #
    # v      : dictionary key value of type float

    # TODO [Low]: Improve descriptions

    if isinstance(v, float):
        return v

    if isinstance(v, int):
        return float(v)

    if isinstance(v, (str, unicode)):
        return float(readStr(v))

    else:
      raise TypeError("key value " + str(v) + " is not an float")



def readDim(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type Dimension

    # TODO [Low]: Improve descriptions

    if isinstance(v, Dimension):
        return v

    if isinstance(v, (str, unicode)):
        dim = readStr(v).strip("[]")
        dims = dim.split()
        return Dimension(*dims)

    else:
      raise TypeError("key value " + str(v) + " has no dimensional type")



def readBool(v):

    # return : interpreted dictionary boolean value
    #
    # v      : dictionary key value of type BoolProxy

    # TODO [Low]: Improve descriptions

    if isinstance(v, BoolProxy):
        return v

    elif isinstance(v, (str, unicode)):
        trueKeys = ["On", "on", "True", "true", "Yes", "yes"]
        falseKeys  = ["Off", "off", "False", "false", "No", "no"]
        for k in trueKeys:
            if readStr(v) == k: return BoolProxy(val=True)
        for k in falseKeys:
            if readStr(v) == k: return BoolProxy(val=False)
        raise TypeError("keyValue " + str(v) + " is no boolean")

    elif isinstance(v, int):
        if v == 1: return BoolProxy(val=True)
        elif v == 0: return BoolProxy(val=False)

    else:
      raise TypeError("key value " + str(v) + " is no boolean")



def readSubDict(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type DictProxy

    # TODO [Low]: Improve descriptions

    if isinstance(v, DictProxy):
        return v

    else:
      raise TypeError("key value " + str(v) + " is no (sub-)dictionary")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

