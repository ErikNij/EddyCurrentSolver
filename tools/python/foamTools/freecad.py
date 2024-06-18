#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module template
# March 2015
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

sys.path.append("/usr/lib/freecad/lib")

import math as m
import numpy as np

import FreeCAD, Sketcher, Draft, Part, PartDesign, Mesh, MeshPart
from FreeCAD import Units, Placement, Matrix, Vector, Rotation
from Part import Line, Circle

from foamTools.math import rotation

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def makeSketch(document, key, orient=None, axis=None, angle=None,
               base=(0.0, 0.0, 0.0), reverse=False):

    name = "Sketch" + key.capitalize()
    label = "sketch_" + key

    if orient:

        if axis: raise KeyError("Axis may not be given if orient is used.")

        if angle: raise KeyError("Angle may not be given if orient is used.")

        if orient == "xy":

            if not reverse: rotation = Rotation(Vector(0.0, 0.0, 1.0), 0)

            else: rotation = Rotation(Vector(-1.0, 0.0, 0.0), 180)

        elif orient == "xz":

            if not reverse: rotation = Rotation(Vector(1.0, 0.0, 0.0), 90)

            else: rotation = Rotation(Vector(0.0, 1.0, 1.0), 180)

        elif orient == "yz":

            if not reverse: rotation = Rotation(Vector(1.0, 1.0, 1.0), 120)

            else: rotation = Rotation(Vector(-1.0, 1.0, 1.0), 240)

        else:

            raise ValueError("Postition " + orient + " is not supported.")

    else:

        rotation = Rotation(Vector(axis), angle)

    sketch = document.addObject("Sketcher::SketchObject", name)
    sketch.Label = label
    sketch.Placement = Placement(Vector(base), rotation)

    return sketch



def sketchCircle(sketchObj, radius, center=(0.0, 0.0, 0.0)):

    document = sketchObj.Document
    document.recompute()

    nConstraints = len(sketchObj.Constraints)
    nGeometry = len(sketchObj.Geometry)

    sketchObj.addGeometry(Circle(Vector(center),
                                 Vector(0.0, 0.0, 1.0), radius))



def sketchPolyLine(sketchObj, v, v2D, fillet=None):

    document = sketchObj.Document
    document.recompute()

    nConstraints = len(sketchObj.Constraints)
    nGeometry = len(sketchObj.Geometry)

    l = len(v) - 1
    L = nGeometry + l

    def V(i): return Vector((v2D[i][0], v2D[i][1], 0.0))

    # Draw lines
    for i in range(l):

        sketchObj.addGeometry(Line(V(v[i]), V(v[i+1])))

    sketchObj.addGeometry(Line(V(v[l]), V(v[0])))

    # Add Constraints (lines need to be drawn first)
    for i in range(l):

        I = nGeometry + i

        sketchObj.addConstraint(
            Sketcher.Constraint('Coincident', I, 2, I+1, 1))

    sketchObj.addConstraint(
        Sketcher.Constraint('Coincident', L, 2, nGeometry, 1))

    if fillet:

        # Add fillets
        for i in range(l+1):

            I = nGeometry + i

            sketchObj.fillet(I, 1, fillet)



def makeFuseBody(key, fuseObjects):

    document = fuseObjects[0].Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    # If only one object is given, fuse with itself
    if len(fuseObjects) == 1: fuseObjects.append(fuseObjects[0])

    body = document.addObject("Part::MultiFuse", name)

    body.Label = label
    body.Shapes = fuseObjects

    return body



def makeCutBody(key, baseObj, toolObj):

    document = baseObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    body = document.addObject("Part::Cut", name)

    body.Label = label
    body.Base = baseObj
    body.Tool = toolObj

    return body



def makeMirrorBody(key, baseObj, normal, base=(0.0, 0.0, 0.0)):

    document = baseObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    body = document.addObject("Part::Mirroring", name)

    body.Label = label
    body.Source = baseObj
    body.Normal = normal
    body.Base = base

    return body



def makeExtrudeBody(key, sketchObj, normal, solid=True, angle=0.0):

    document = sketchObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    if isinstance(normal, (int, float)):

        rotAxis = np.array(sketchObj.Placement.Rotation.Axis)
        rotTheta = sketchObj.Placement.Rotation.Angle

        direction = tuple(rotation(rotAxis, rotTheta,
                                   np.array([0.0, 0.0, normal])))

    else:

        direction = tuple(np.array(normal).copy())

    body = document.addObject("Part::Extrusion", name)

    body.Label = label
    body.Base = sketchObj
    body.Dir = Vector(direction)
    body.Solid = solid
    body.TaperAngle = angle

    return body



def makeDoubleExtrudeBody(key, sketchObj, normal, solid=True, angle=0.0):

    document = sketchObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    if isinstance(normal, (int, float)):

        rotAxis = np.array(sketchObj.Placement.Rotation.Axis)
        rotTheta = sketchObj.Placement.Rotation.Angle

        direction = tuple(rotation(rotAxis, rotTheta,
                                   np.array([0.0, 0.0, normal/2.0])))

    else:

        direction = tuple(np.array(normal).copy()/2.0)

    front = makeExtrudeBody(key + "_front", sketchObj,
                              normal=direction, solid=solid, angle=angle)

    direction = tuple(-np.array(direction))

    back = makeExtrudeBody(key + "_back", sketchObj,
                              normal=direction, solid=solid, angle=angle)

    body = makeFuseBody(key, [front, back])

    return body



def makeRevolveBody(key, sketchObj, angle=360.0, axis=(0.0, 0.0, 1.0),
                    base=(0.0, 0.0, 0.0), solid=True):

    document = sketchObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    body = document.addObject("Part::Revolution", name)

    body.Label = label
    body.Source = sketchObj
    body.Axis = Vector(axis)
    body.Base = Vector(base)
    body.Angle = angle
    body.Solid = solid

    return body



def makeDoubleRevolveBody(key, sketchObj, angle=360.0, axis=(0.0, 0.0, 1.0),
                          base=(0.0, 0.0, 0.0), solid=True):

    document = sketchObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    front = makeRevolveBody(key + "_front", sketchObj, angle=angle/2.0,
                              axis=axis, base=base, solid=solid)

    back = makeRevolveBody(key + "_back", sketchObj, angle=-angle/2.0,
                              axis=axis, base=base, solid=solid)

    body = makeFuseBody(key, [front, back])

    return body



def makeOrthoArrayBody(key, baseObj, xvector, xnum,
                       yvector=(0.0, 0.0, 0.0), ynum=1,
                       zvector=(0.0, 0.0, 0.0), znum=1, fuse=False):

    document = baseObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    body = Draft.makeArray(baseObj, Vector(xvector), Vector(yvector),
                           xnum, ynum, name=name)

    body.Label = label
    body.ArrayType = str("ortho")
    body.Base = baseObj
    body.IntervalX = Vector(xvector)
    body.IntervalY = Vector(yvector)
    body.IntervalZ = Vector(zvector)
    body.NumberX = xnum
    body.NumberY = ynum
    body.NumberZ = znum
    body.Fuse = fuse


    document.recompute()

    return body



def makePolarArrayBody(key, baseObj, totalnum, center=(0.0, 0.0, 0.0),
                       axis=(0.0, 0.0, 1.0), iaxis=(0.0, 0.0, 0.0),
                       totalangle=360.0, fuse=False):

    document = baseObj.Document
    document.recompute()

    name = "Body" + key.capitalize()
    label = "body_" + key

    document.recompute()

    body = Draft.makeArray(baseObj, center, totalangle, totalnum, name=name)

    body.Label = label
    body.ArrayType = str("polar")
    body.Base = baseObj
    body.Center = center
    body.Axis = axis
    body.IntervalAxis = iaxis
    body.Angle = totalangle
    body.NumberPolar = totalnum
    body.Fuse = fuse

    return body



def faceShell(shd):

    dList = list([shd]) if isinstance(shd, tuple) else shd

    if not isinstance(dList, list):

        raise ValueError("Shell data must be list of tuples.")

    fList = list()

    for i in range(len(dList)):

        d = dList[i]

        if not isinstance(d, tuple):

            raise ValueError("Shell dictionary data must be tuple.")

        body = d[0]
        labels = d[1]

        if i == 0:

            document = body.Document
            document.recompute()

        lList = [labels] if isinstance(labels, int) else labels

        if not isinstance(lList, list):

            raise ValueError("Face labels must be list.")

        for l in lList:

            if not isinstance(l, int):

                raise ValueError("Face label must be integer.")

            eval("fList.append(body.Shape.Face" + str(l) + ")")

    return Part.Shell(fList)



def makeFaceShell(document, key, shd):

    name = "Shell" + key.capitalize()
    label = "shell_" + key

    shell = document.addObject("Part::Feature", name)

    shell.Label = label
    shell.Shape = faceShell(shd)

    return shell



def exportMeshes(objects, dir, prefix, tol=0.1, scale=1.0):

    document = objects[0].Document
    document.recompute()

    S = Matrix()
    S.scale(scale, scale, scale)

    for i in range(len(objects)):

        o = objects[i]

        mesh = Mesh.Mesh(o.Shape.tessellate(tol))
        #mesh = MeshPart.meshFromShape(Shape=o.Shape, MaxLength=10)

        mesh.transform(S)

        name = prefix + "_" + o.Label + ".stl"

        mesh.write(dir + "/" + name)

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    pass

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

