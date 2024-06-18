#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# June 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_d = 0.0032
geo_r = 0.00736 + geo_d/2.0
geo_z = 0.00343

top_n = 2
top_z = 0.00465

bottom_n = 3
bottom_z = 0.00465

arc_n = 36
arc_alpha = (2.0*m.pi)/arc_n

bundle_n = 36
bundle_alpha = (2.0*m.pi)/bundle_n

# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

# --------------------------------------------------------------------------- #

for top in range(top_n):

    # File name base
    fileNameBase = 'coil_top'

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    for bundle in range(bundle_n):

        # Loop
        for n in range(arc_n):

            br = geo_d/2.0 * m.cos(bundle*bundle_alpha)
            bz = geo_d/2.0 * m.sin(bundle*bundle_alpha)

            x =  (geo_r + br) * m.cos(n*arc_alpha)
            y =  (geo_r + br) * m.sin(n*arc_alpha)
            z =  geo_z + top*top_z + bz

            points.append([x, y, z])

        shift = bundle*arc_n

        # Points string
        for pointi in range(len(points)-shift):

            pointI = pointi + shift

            point = points[pointI]

            p += i(1, '('+str(point[0])+' '+str(point[1])+' '+str(point[2])+')')

        # Edges string
        for pointi in range(len(points)-1-shift):

            pointI = pointi + shift

            e += i(1, '('+str(pointI)+' '+str(pointI+1)+')')

        e += i(1, '('+str(len(points)-1)+' '+str(shift)+')')

    # Write
    fileName = fileNameBase + str(top) + '.eMesh'
    with open(fileGetPath(fileName),'w') as f:

        f.write(objectHeader(fileName, 'featureEdgeMesh'))

        f.write(i(0, '// points:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(p)
        f.write(i(0, ')\n'))

        f.write(i(0, '// edges:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(e)
        f.write(i(0, ')\n'))

        f.write(objectFooter())


for bottom in range(bottom_n):

    # File name base
    fileNameBase = 'coil_bottom'

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    for bundle in range(bundle_n):

        # Loop
        for n in range(arc_n):

            br = geo_d/2.0 * m.cos(bundle*bundle_alpha)
            bz = geo_d/2.0 * m.sin(bundle*bundle_alpha)

            x =  (geo_r + br) * m.cos(n*arc_alpha)
            y =  (geo_r + br) * m.sin(n*arc_alpha)
            z =  -geo_z - bottom*bottom_z + bz

            points.append([x, y, z])

        shift = bundle*arc_n

        # Points string
        for pointi in range(len(points)-shift):

            pointI = pointi + shift

            point = points[pointI]

            p += i(1, '('+str(point[0])+' '+str(point[1])+' '+str(point[2])+')')

        # Edges string
        for pointi in range(len(points)-1-shift):

            pointI = pointi + shift

            e += i(1, '('+str(pointI)+' '+str(pointI+1)+')')

        e += i(1, '('+str(len(points)-1)+' '+str(shift)+')')

    # Write
    fileName = fileNameBase + str(bottom) + '.eMesh'
    with open(fileGetPath(fileName),'w') as f:

        f.write(objectHeader(fileName, 'featureEdgeMesh'))

        f.write(i(0, '// points:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(p)
        f.write(i(0, ')\n'))

        f.write(i(0, '// edges:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(e)
        f.write(i(0, ')\n'))

        f.write(objectFooter())


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
