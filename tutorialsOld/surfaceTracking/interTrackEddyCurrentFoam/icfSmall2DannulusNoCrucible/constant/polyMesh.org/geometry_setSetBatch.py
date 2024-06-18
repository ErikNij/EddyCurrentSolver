#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.blockMeshDict import ioBase
import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

from geometry_parameters import *

def s(val):

    return str(geo_scale * val)

# --------------------------------------------------------------------------- #
# --- setSetBatch ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

ssb = ioBase("setSetBatch")

ssb.write("### Geometry")
ssb.line()

# --------------------------------------------------------------------------- #

ssb.write("# All")
ssb.line()

ssb.write("cellSet cellSet_geometry_all                 new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + str(-1) + " " + str(-1) + " " + str(-1) + " )", end=" ")
ssb.write("( " + str(1) + " " + str(1) + " " + str(1) + " )")
ssb.line()



ssb.write("# Fluid")
ssb.line()

ssb.write("cellSet cellSet_geometry_fluid_dynamic       new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(0) + " " + str(-1) + " " + s(0) + " )", end=" ")
ssb.write("( " + s(geo_x0) + " " + str(1) + " " + s(geo_z2) + " )")
ssb.line()

ssb.write("cellSet cellSet_geometry_fluid               new", end=" ")
ssb.write("setToCell cellSet_geometry_fluid_dynamic")
ssb.line()



ssb.write("# Crucible")
ssb.line()

ssb.write("cellSet cellSet_geometry_crucible_dynamic    new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(geo_x0) + " " + str(-1) + " " + s(0) + " )", end=" ")
ssb.write("( " + s(geo_x1) + " " + str(1) + " " + s(geo_z3) + " )")
ssb.line()

ssb.write("cellSet cellSet_geometry_crucible_static     new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(0) + " " + str(-1) + " " + s(geo_x0-geo_x1) + " )", end=" ")
ssb.write("( " + s(geo_x1) + " " + str(1) + " " + s(geo_z3) + " )")
ssb.write("#cellSet cellSet_geometry_crucible_static     add", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(geo_x0) + " " + str(-1) + " " + s(geo_z3) + " )", end=" ")
ssb.write("( " + s(geo_x1) + " " + str(1) + " " + s(geo_z3) + " )")
ssb.write("#cellSet cellSet_geometry_crucible_static_old remove")
ssb.line()

ssb.write("cellSet cellSet_geometry_crucible            new", end=" ")
ssb.write("setToCell cellSet_geometry_crucible_dynamic")
ssb.write("cellSet cellSet_geometry_crucible            add", end=" ")
ssb.write("setToCell cellSet_geometry_crucible_static")
ssb.write("cellSet cellSet_geometry_crucible_old        remove")
ssb.line()



ssb.write("# Gap")
ssb.line()

ssb.write("cellSet cellSet_geometry_gap_dynamic         new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(geo_x1) + " " + str(-1) + " " + s(0) + " )", end=" ")
ssb.write("( " + s(geo_x2) + " " + str(1) + " " + s(geo_z3) + " )")
ssb.line()

ssb.write("#cellSet cellSet_geometry_gap_static          new", end=" ")
ssb.write("boxToCell", end=" ")
ssb.write("( " + s(geo_x1) + " " + str(-1) + " " + s(geo_z3) + " )", end=" ")
ssb.write("( " + s(geo_x2) + " " + str(1) + " " + s(geo_z3) + " )")
ssb.line()

ssb.write("cellSet cellSet_geometry_gap                 new", end=" ")
ssb.write("setToCell cellSet_geometry_gap_dynamic")
ssb.write("cellSet cellSet_geometry_gap                 add", end=" ")
ssb.write("setToCell cellSet_geometry_gap_static")
ssb.write("cellSet cellSet_geometry_gap_old             remove")
ssb.line()



#cellSet cellSet_geometry_crucibleDynamic  new boxToCell ( 0.08500 -1.00000  0.00000) ( 0.08600  1.00000  0.38300)

#cellSet cellSet_geometry_crucibleStatic   new boxToCell ( 0.00000 -1.00000 -0.00100) ( 0.08600  1.00000  0.00000)
#cellSet cellSet_geometry_crucibleStatic   add boxToCell ( 0.08500 -1.00000  0.38300) ( 0.08600  1.00000  0.38300)

#cellSet cellSet_geometry_crucible         new setToCell cellSet_geometry_crucibleDynamic
#cellSet cellSet_geometry_crucible         add setToCell cellSet_geometry_crucibleStatic

# --------------------------------------------------------------------------- #

ssb.write("### End")
ssb.line()

ssb.write("quit")
ssb.line()


# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
